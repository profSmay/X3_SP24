#region imports
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

from copy import deepcopy as dc
#endregion

#region class definitions
class StateDataForPlotting:
    """
    I'm making this class for easy storage of data for plotting.
    """
    def __init__(self):
        self.T = []
        self.P = []
        self.h = []
        self.u = []
        self.s = []
        self.v = []

    def clear(self):
        self.T.clear()
        self.P.clear()
        self.h.clear()
        self.u.clear()
        self.s.clear()
        self.v.clear()

    def add(self, vals):
        T, P, u, h, s, v = vals
        self.T.append(T)
        self.P.append(P)
        self.h.append(h)
        self.u.append(u)
        self.s.append(s)
        self.v.append(v)

    def getAxisLabel(self, W='T', Units=None):
        Units = Units if Units is not None else units()
        w=W.lower()
        if w == 't':
            return Units.TPlotUnits
        if w == 'h':
            return Units.hPlotUnits
        if w == 'u':
            return Units.uPlotUnits
        if w == 's':
            return Units.sPlotUnits
        if w == 'v':
            return Units.vPlotUnits
        if w == 'p':
            return Units.PPlotUnits

    def getDataCol(self, W='T'):
        w=W.lower()
        if w=='t':
            return self.T
        if w=='h':
            return self.h
        if w=='u':
            return self.u
        if w=='s':
            return self.s
        if w=='v':
            return self.v
        if w=='p':
            return self.P

class stateProps():
    """
    for storage and retrieval of a thermodynamic state
    T, P, u, h, s, v
    """
    def __init__(self):
        self.name = None
        self.T = None
        self.P = None
        self.h = None
        self.u = None
        self.s = None
        self.v = None

    # this is overloading the multiply operator.  Allows me to multiply a scalar or do a dot product (i.e., b=s*a or c=b*a)
    def __mul__(self, other):
        if type(other) in (float, int):
            b=stateProps()
            b.h*=other
            b.u*=other
            b.s*=other
            b.v*=other
            return b

    # this is overloading the __rmul__ operator so that s*Pt works.
    def __rmul__(self,other):
        return self*other

    # this is overloading the division operator.  Allows me to divide by a scalar (i.e., b=a/s)
    def __truediv__(self, other):
        if type(other) in (float, int):
            b = stateProps()
            b.h /= other
            b.u /= other
            b.s /= other
            b.v /= other
            return b

    def ConvertStateData(self, SI=True, mass=False, total=False, n=1.0, MW=1.0, Units=None):
        UC=Units if Units is not None else units()
        UC.set(SI=SI, mass=mass, total=total)
        mCF = 1.0 if SI else UC.CF_Mass
        TCF = 1.0 if SI else UC.CF_T
        PCF = 1.0 if SI else UC.CF_P
        vCF = 1.0 if SI else UC.CF_v  # convert m^3/mol to ft^3/lbmol
        uCF = 1.0 if SI else UC.CF_e  # cpmvert J/mol to Btu/lbmol
        hCF = 1.0 if SI else UC.CF_e
        sCF = 1.0 if SI else UC.CF_s
        nCF = 1.0 if SI else UC.CF_n  # convert mol to lbmol
        if mass:
            mCF/=MW
            vCF/=MW
            uCF/=MW
            hCF/=MW
            sCF/=MW
        elif total:
            vCF*=n*nCF
            uCF*=n*nCF
            hCF*=n*nCF
            sCF*=n*nCF

        self.P*=PCF
        self.T*=TCF
        self.h*=hCF
        self.u*=uCF
        self.v*=vCF
        self.s*=sCF

    def getVal(self, name='T'):
        n=name.lower()
        if n == 't':
            return self.T
        if n == 'h':
            return self.h
        if n == 'u':
            return self.u
        if n == 's':
            return self.s
        if n == 'v':
            return self.v
        if n == 'p':
            return self.P

    def print(self, units=None):
        U = units if units is not None else units()
        if self.name is not None:
            print(self.name)
        print('v={:0.4f} {}.'.format(self.v, self.U.vUnits))
        print('u={:0.4f} {}'.format(self.u,self.U.uUnits))
        print('h={:0.4f} {}'.format(self.h, self.U.hUnits))
        print('s={:0.4f} {}'.format(self.s, self.U.sUnits))

class units():
    """
    For air, I'm assuming the default units are on a molar basis.
    """
    def __init__(self):
        self.SI=True
        # default set of units
        self.sUnits='J/mol*k'
        self.uUnits='J/mol'
        self.vUnits='m^3/mol'
        self.VUnits='m^3'
        self.hUnits=self.uUnits
        self.mUnits='kg'
        self.TUnits='K'
        self.PUnits='Pa'
        self.EUnits='J'
        
        # conversion factors
        self.CF_E = 1.0/1055.06  # J to Btu
        self.CF_Length = 3.28084  # m to ft
        self.CF_V = self.CF_Length**3.0  # m^3 to ft^3
        self.CF_P = 1.0/101325  # Pa to atm
        self.CF_Mass = 2.20462  # kg to lb
        self.CF_T = 9.0/5.0  # K to R
        self.CF_n = 1/453.59 # mol to lbmol
        self.CF_v = self.CF_V/self.CF_n  # m^3/mol to ft^3/lbmol
        self.CF_e = self.CF_E/self.CF_n  # J/mol to Btu/lbmol
        self.CF_s = self.CF_e/(self.CF_n*self.CF_T)  #J/mol*K to Btu/lbmol*R

        self.setPlotUnits()

    def set(self, SI=True, mass=False, total=False):
        self.changed = not self.SI == SI
        self.SI=SI
        if SI:
            self.sUnits = 'J/{}k'.format('' if total else ('kg*' if mass else 'mol*'))
            self.uUnits = 'J{}'.format('' if total else ('/kg' if mass else '/mol'))
            self.vUnits = 'm^3{}'.format('' if total else ('/kg' if mass else '/mol'))
            self.VUnits = 'm^3'
            self.hUnits = self.uUnits
            self.mUnits = 'kg'
            self.TUnits = 'K'
            self.PUnits = 'Pa'
            self.EUnits = 'J'
        else:
            self.sUnits = 'BTU/{}R'.format('' if total else ('lb*' if mass else 'lbmol*'))
            self.uUnits = 'BTU{}'.format('' if total else ('/lb' if mass else '/lbmol'))
            self.vUnits = 'ft^3{}'.format('' if total else ('/lb' if mass else '/lbmol'))
            self.VUnits = 'ft^3'
            self.hUnits = self.uUnits
            self.mUnits = 'lb'
            self.TUnits = 'R'
            self.PUnits = 'Atm'
            self.EUnits = 'Btu'

        self.setPlotUnits(SI=SI, mass=mass, total=total)
    
    def setPlotUnits(self, SI=True, mass=True, total=False):
        if SI:
            self.PPlotUnits = r'P $\left(Pa\right)$'
            self.TPlotUnits = r'T $\left(K\right)$'
            if total:
                self.sPlotUnits = r'S $\left(\frac{J}{K}\right)$'
                self.uPlotUnits = r'U $\left(J\right)$'
                self.hPlotUnits = r'H $\left(J\right)$'
                self.vPlotUnits = r'V $\left(m^3\right)$'
            elif mass:
                self.sPlotUnits = r's $\left(\frac{J}{kg*K}\right)$'
                self.uPlotUnits = r'u $\left(\frac{J}{kg}\right)$'
                self.hPlotUnits = r'h $\left(\frac{J}{kg}\right)$'
                self.vPlotUnits = r'v $\left(\frac{m^3}{kg}\right)$'
            else:
                self.sPlotUnits = r'$\bar{s} \left(\frac{J}{mol*K}\right)$'
                self.uPlotUnits = r'$\bar{u} \left(\frac{J}{mol}\right)$'
                self.hPlotUnits = r'$\bar{h} \left(\frac{J}{mol}\right)$'
                self.vPlotUnits = r'$\bar{v} \left(\frac{m^3}{mol}\right)$'
        else:
            self.PPlotUnits = r'P $\left(atm\right)$'
            self.TPlotUnits = r'T $\left(^{o}R\right)$'
            if total:
                self.sPlotUnits = r'S $\left(\frac{Btu}{^{o}R}\right)$'
                self.uPlotUnits = r'U $\left(Btu\right)$'
                self.hPlotUnits = r'H $\left(Btu\right)$'
                self.vPlotUnits = r'V $\left(ft^3\right)$'
            elif mass:
                self.sPlotUnits = r's $\left(\frac{Btu}{lb\cdot^{o}R}\right)$'
                self.uPlotUnits = r'u $\left(\frac{Btu}{lb}\right)$'
                self.hPlotUnits = r'h $\left(\frac{Btu}{lb}\right)$'
                self.vPlotUnits = r'v $\left(\frac{ft^3}{lb}\right)$'
            else:
                self.sPlotUnits = r'$\bar{s} \left(\frac{Btu}{lb_{mol}\cdot^{o}R}\right)$'
                self.uPlotUnits = r'$\bar{u} \left(\frac{Btu}{lb_{mol}}\right)$'
                self.hPlotUnits = r'$\bar{h} \left(\frac{Btu}{lb_{mol}}\right)$'
                self.vPlotUnits = r'$\bar{v} \left(\frac{ft^3}{lb_{mol}}\right)$'

    #region temperature conversion formulas
    def T_RtoK(self, T):
        return T*5.0/9.0

    def T_FtoC(self, T):
        return (T-32.0)*5.0/9.0

    def T_RtoF(self, T):
        return T-459.67

    def T_FtoK(self, T):
        return self.T_RtoK(self.T_FtoR(T))

    def T_CtoK(self, T):
        return T+273.15

    def T_CtoF(self, T):
        return T * 9.0 / 5.0 + 32

    def T_KtoC(self, T):
        return T - 273.15

    def T_KtoR(self, T):
        return T*9/5

    def T_FtoR(self, T):
        return T + 459.67
    #endregion

class air():
    def __init__(self):
        """
        Air as an ideal gas.
        I choose to always specify air in molar metric units.
        The standard state is T=0C and P=1 atm (101.325kPa)
        So u0=0, h0=0, s0=0, v0=RT/P
        :param P: Pressure in Pa
        :type P: float
        :param T: Temperature in K
        :type T: float
        :param V: Volume in m^3/mol
        :type V: float
        """
        self.RBar = 8.3145 # J/mol*K or kJ/kmol*K
        self.MW = 28.97 # kg/kmol or g/mol or lb/lbmol
        self.R=self.RBar/self.MW  # kJ/kg*K or J/g*K
        #region set standard state properties
        self.StandardState=stateProps()
        self.StandardState.P = 101325.0 # P in Pa
        self.StandardState.T = 273.15 # T in K
        self.StandardState.v = self.RBar*self.StandardState.T/self.StandardState.P # v in m^3/mol
        self.StandardState.u=0
        self.StandardState.h=0
        self.StandardState.s=0
        #endregion
        self.State=stateProps()
        self.n = 1.0  # moles
        self.m=self.n*self.MW/1000.0  # mass in kg

    def cv(self, T):
        return self.cp(T)-self.RBar

    def cp(self, T):
        """
        For air as an ideal gas, cp is a function of temperature as given by:
        cp=Rbar(a+b*T+c*T**2+d*T**3+e*T**4)
        :param T: is Temperature in K
        :type T: float
        :return: molar specific heat in units of kJ/kg
        :rtype: float
        """
        TLowRange=1630.0
        a = 3.653 if T<TLowRange else 2.753
        b = -1.337E-3 if T<TLowRange else 0.002
        c = 3.294E-6 if T<TLowRange else -1.0E-6
        d = -1.913E-9  if T<TLowRange else 3.0E-10
        e = 0.2763E-12  if T<TLowRange else -3.0E-14
        return self.RBar*(a+b*T+c*T**2+d*T**3+e*T**4)

    def deltau(self, T1=None, T2=None):
        """
        To calculate changes in molar internal energy for air as an ideal gas u=u(T)
        cv=du/dT|v -> delta u=int((cv)dT, T1, T2)
        :param T1: Temperature 1 in K
        :type T1: float
        :param T2: Temperature 2 in K
        :type T2: float
        :return: deltau in kJ/kmol or J/mol
        :rtype: float
        """
        if T1 is None:
            T1=self.StandardState.T
        if T2 is None:
            T2=self.StandardState.T
        return quad(self.cv,T1,T2)[0]

    def deltah(self, T1=None, T2=None):
        """
        To calculate changes in molar internal energy for air as an ideal gas u=u(T)
        cp=dh/dT|p -> delta h=int((cp)dT, T1, T2)
        :param T1: temperature 1 in K
        :type T1: float
        :param T2: temperature 2 in K
        :type T2: float
        :return: delta h in kJ/kmol or J/mol
        :rtype: float
        """
        if T1 is None:
            T1=self.StandardState.T
        if T2 is None:
            T2 = self.StandardState.T
        return quad(self.cp,T1,T2)[0]

    def deltas_tv(self, T1=None, T2=None, V1=None, V2=None):
        """
        For calculating changes in molar entropy for air as an ideal gas s=s(T,V)
        Tds=du+Pdv -> delta s = int(cv/T*dT, T1, T2)+R ln(V2/V1)
        :param T1: Temperature 1 in K
        :type T1: float
        :param T2:  Temperature 2 in K
        :type T2: float
        :param V1:  Volume 1
        :type V1: float
        :param V2:  Volume 2
        :type V2: float
        :return: delta s in J/mol*K
        :rtype: float
        """
        if T1 is None:
            T1=self.StandardState.T
        if T2 is None:
            T2 = self.StandardState.T
        if V1 is None:
            V1 = self.StandardState.v
        if V2 is None:
            V2 = self.StandardState.v
        fn=lambda T: 0 if T==0 else self.cv(T)/T
        deltaS=quad(fn,T1,T2)[0]
        deltaS+=self.RBar*math.log(V2/V1)
        return deltaS

    def deltas_tp(self, T1=None, T2=None, P1=None, P2=None):
        """
        For calculating changes in molar entropy for air as an ideal gas s=s(T,V)
        Tds=dh-vdP -> delta s = int(cp/T*dT, T1, T2)-R ln(P2/P1)
        :param T1: Temperature 1 in K
        :type T1: float
        :param T2:  Temperature 2 in K
        :type T2: float
        :param P1:  pressure 1 in Pa
        :type P1: float
        :param P2:  pressure 2 in Pa
        :type P2: float
        :return: delta s in J/mol*K
        :rtype: float
        """
        if T1 is None:
            T1 = self.StandardState.T
        if T2 is None:
            T2 = self.StandardState.T
        if P1 is None:
            P1 = self.StandardState.P
        if P2 is None:
            P2 = self.StandardState.P

        fn=lambda T: 0 if T==0.0 else self.cp(T)/T
        deltaS=quad(fn,T1,T2)[0]
        deltaS+=self.RBar*math.log(P1/P2)
        return deltaS

    def set(self, P=None, T=None, v=None, h=None, u=None, s=None, name=None):
        """
        This allows me to set two properties and calculate the state of the air
        :param pressure: in Pa
        :param T: Temperature in K
        :param v: specific volume in m^3/mol
        :param u: specific internal energy in J/mol
        :param h: specific enthalpy in J/mol
        :param s: specific entropy in J/mol*K
        :param name: a convenient name
        :return: a deep copy of the calculated state
        """
        self.State.P = P  # pressure - Pa
        self.State.T = T  # Temperature - K
        self.State.v = v  # specific volume - m^3/mol
        self.State.h = h  # specific enthalpy - J/mol
        self.State.u = u  # specific internal energy - J/mol
        self.State.s = s  # entropy - J/(mol*K)
        self.State.name=name
        if T == None and P==None and u==None and v == None and h == None and s == None:
            return
        else:
            self.calc()
        return dc(self.State)  # need to deep copy so not passing just a reference back

    def calc(self):
        '''
        To calculate the state of ideal gas air, we use the ideal gas law and specific heat functions relative to
        the standard state of T=0C, P=101.325 kPa where u=0, h=0, s=0, v=vo by declaration
        In the general case, we have 6 thermodynamic properties (dof), but 2 are specified:  6!/4!2!=15 permutations.
        P: T, u, v, h, s
        T: u, v, h, s  (because u & h are only dependent on T for an ideal gas, specifying T+u or T+h does not work)
        u: v, h, s  (because u & h are only dependent on T for an ideal gas, specifying u+h does not work)
        v: h, s
        h: s
        :return: a deep copy of the state in specific molar properties
        '''
        # 1. need to determine which two properties are known
        # 2. calculate all the other thermodynamic properties
        State=stateProps()
        #region case 1. P,T
        if self.State.P is not None and self.State.T is not None:
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.u=self.deltau(T2=self.State.T)
            self.State.h=self.deltah(T2=self.State.T)
            self.State.s=self.deltas_tp(T2=self.State.T, P2=self.State.P)
        #endregion
        #region case 2. P,u
        elif self.State.P is not None and self.State.u is not None:
            fn = lambda T: self.deltau(T2=T[0])-self.State.u
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.h=self.deltah(T2 = self.State.T)
            self.State.s=self.deltas_tp(T2=self.State.T,P2=self.State.P)
        #endregion
        #region case 3. P,v
        elif self.State.P is not None and self.State.v is not None:
            self.State.T=self.State.v*self.State.P/self.RBar
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.u=self.deltau(T2=self.State.T)
            self.State.h=self.deltah(T2=self.State.T)
            self.State.s=self.deltas_tp(T2=self.State.T,P2=self.State.P)
        #endregion
        #region case 4. P,h
        elif self.State.P is not None and self.State.h is not None:
            fn = lambda T: self.deltah(T2=T[0])-self.State.h
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.u=self.deltau(T2=self.State.T)
            self.State.s=self.deltas_tp(T2=self.State.T,P2=self.State.P)
        #endregion
        #region case 5. P,s
        elif self.State.P is not None and self.State.s is not None:
            fn = lambda T: self.deltas_tp(T2=T[0], P2=self.State.P)-self.State.s
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.u=self.deltau(T2=self.State.T)
            self.State.h=self.deltah(T2=self.State.T)
        #endregion
        #region case 6. T,u  # T & u not independent
        #endregion
        #region case 7. T,v
        elif self.State.T is not None and self.State.v is not None:
            self.State.P=self.State.T*self.RBar/self.State.v
            self.State.u=self.deltau(T2=self.State.T)
            self.State.h=self.deltah(T2=self.State.T)
            self.State.s=self.deltas_tp(T2=self.State.T,P2=self.State.P)
        #endregion
        #region case 8. T,h # T & h not independent
        #endregion
        #region case 9. T,s
        elif self.State.T is not None and self.State.s is not None:
            fn = lambda P: self.deltas_tp(T2=self.State.T, P2=P[0])-self.State.s
            r = fsolve(fn, np.array([50]))
            self.State.P=r[0]
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.u=self.deltau(T2=self.State.T)
            self.State.h=self.deltah(T2=self.State.T)
        #endregion
        #region case 10. T,v
        elif self.State.u is not None and self.State.v is not None:
            fn = lambda T: self.deltau(T2=T[0])-self.State.u
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            self.State.P=self.State.T*self.RBar/self.State.v
            self.State.h=self.deltah(T2=self.State.T)
            self.State.s=self.deltas_tp(T2=self.State.T,P2=self.State.P)
        #endregion
        #region case 11. u,h # u & h not independent
        #endregion
        #region case 12. T,s
        elif self.State.u is not None and self.State.s is not None:
            fn = lambda T: self.deltau(T2=T[0])-self.State.u
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            fn = lambda P: self.deltas_tp(T2=self.State.T, P2=P)-self.State.s
            r = fsolve(fn, np.array([50]))
            self.State.P=r[0]
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.h=self.deltah(T2=self.State.T)
        #endregion
        #region case 13. v,h
        elif self.State.v is not None and self.State.h is not None:
            fn = lambda T: self.deltah(T2=T[0])-self.State.h
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            self.State.P=self.State.T*self.RBar/self.State.v
            self.State.u=self.deltau(T2=self.State.T)
            self.State.s=self.deltas_tp(T2=self.State.T, P2=self.State.P)
        #endregion
        #region case 14. v,s
        elif self.State.v is not None and self.State.s is not None:
            fn = lambda T: self.deltas_tv(T2=T[0], V2=self.State.v)-self.State.s
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            self.State.P = self.RBar * self.State.T / self.State.v
            self.State.h = self.deltah(T2=self.State.T)
            self.State.u = self.deltau(T2=self.State.T)
        #endregion
        #region case 15. h,s
        elif self.State.h is not None and self.State.s is not None:
            fn = lambda T: self.deltah(T2=T[0])-self.State.h
            r = fsolve(fn, np.array([50]))
            self.State.T=r[0]
            fn = lambda P: self.deltas_tp(T2=self.State.T, P2=P)-self.State.s
            r = fsolve(fn, np.array([50]))
            self.State.P=r[0]
            self.State.v=self.RBar*self.State.T/self.State.P
            self.State.u=self.deltau(T2=self.State.T)
        #endregion

    def getSummary_MassBasis(self, units=None):
        UC=units if units is not None else units()
        mCF=1.0 if UC.SI else UC.CF_Mass
        TCF=1.0 if UC.SI else UC.CF_T
        PCF=1.0 if UC.SI else UC.CF_P
        vCF=1.0 if UC.SI else UC.CF_V
        uCF=1.0 if UC.SI else UC.CF_E
        hCF=1.0 if UC.SI else UC.CF_E
        sCF=1.0 if UC.SI else UC.CF_S

        stTmp=''
        stTmp+='T={:0.2f} {}\n'.format(self.State.T*TCF, UC.TUnits)
        stTmp+='P={:0.3f} {}\n'.format(self.State.P*PCF/1000.0, UC.PUnits)
        stTmp+='v={:0.4f} {}\n'.format(self.State.v*vCF*1000.0/self.MW, UC.vUnits)
        stTmp+='u={:0.4f} {}\n'.format(self.State.u*uCF/self.MW, UC.uUnits)
        stTmp+='h={:0.4f} {}\n'.format(self.State.h*hCF/self.MW, UC.hUnits)
        stTmp+='s={:0.4f} {}'.format(self.State.s*sCF/self.MW, UC.sUnits)
        return stTmp
        
    def print_MassBasis(self):
        print(self.getSummary_MassBasis())

    def getSummary_Extensive(self, units=None):
        UC=units if units is not None else units()
        mCF = 1.0 if UC.SI else UC.CF_Mass
        TCF = 1.0 if UC.SI else UC.CF_T
        PCF = 1.0 if UC.SI else UC.CF_P
        vCF = 1.0 if UC.SI else UC.CF_V
        uCF = 1.0 if UC.SI else UC.CF_E
        hCF = 1.0 if UC.SI else UC.CF_E
        sCF = 1.0 if UC.SI else UC.CF_S

        stTmp = ''
        stTmp += 'T={:0.2f} {}\n'.format(self.n*self.State.T * TCF, UC.TUnits)
        stTmp += 'P={:0.3f} {}\n'.format(self.n*self.State.P * PCF / 1000.0, UC.PUnits)
        stTmp += 'v={:0.4f} {}\n'.format(self.n*self.State.v * vCF * 1000.0 , UC.vUnits)
        stTmp += 'u={:0.4f} {}\n'.format(self.n*self.State.u * uCF , UC.uUnits)
        stTmp += 'h={:0.4f} {}\n'.format(self.n*self.State.h * hCF , UC.hUnits)
        stTmp += 's={:0.4f} {}'.format(self.n*self.State.s * sCF, UC.sUnits)
        return stTmp

    def print_Extensive(self):
        ext=self.State*self.n
        print('T={:0.2f} {}'.format(ext.T, 'K'))
        print('P={:0.3f} {}'.format(ext.P/1000.0, 'kPa'))
        print('v={:0.4f} {}.'.format(ext.v, 'm^3'))
        print('u={:0.4f} {}'.format(ext.u, 'kJ'))
        print('h={:0.4f} {}'.format(ext.h, 'kJ'))
        print('s={:0.4f} {}'.format(ext.s, 'kJ/K'))
#endregion


def main():
    a=air()
    a.set(P=a.StandardState.P, T=200)
    a.print_Extensive()

if __name__ == "__main__":
    main()