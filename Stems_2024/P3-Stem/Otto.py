#region imports
from Air import *
from matplotlib import pyplot as plt
from PyQt5 import QtWidgets as qtw
import sys
#endregion

#region class definitions
class ottoCycleModel():
    def __init__(self, p_initial=1000.0, v_cylinder=1.0, t_initial=298, t_high=1500.0, ratio=6.0, name='Air Standard Otto Cycle'):
        """
        Constructor for an air standard otto cycle.  The Otto has 4 primary states and consists of four thermodynamic
        processes:
        1. Isentropic compression from: v1, T1, P1 to v2, T2, P2 (Note v2=v1/C.R.)
        2. Constant volume heat addition:  v3=v2
        3. Isentropic expansion (power stroke): v3=v1
        4. Constant volume heat rejection.
        Compression stroke work = (u2-u1)
        Power stroke work = (u3-u4)
        Heat in = (u3-u2)
        Heat out = (u4-u1)
        :param p_initial: Pressure in Pa
        :type p_initial: float
        :param v_cylinder: Volume in m^3
        :type v_cylinder: float
        :param t_initial: Initial Temperature in K
        :type t_initial: float
        :param t_high: High Temperature in K
        :type t_high: float
        :param ratio: Compression ratio
        :type ratio: float
        :param name: a name
        :type name: string
        """
        self.units=units()
        self.units.SI=False
        self.air = air()  # the working fluid
        self.air.set(P=p_initial, T=t_initial)  # initial state if fixed at p_initial, t_initial
        self.p_initial=p_initial
        self.T_initial=t_initial
        self.T_high=t_high
        self.Ratio=ratio  # the compression ratio V_BDC/V_TDC
        self.V_Cylinder=v_cylinder
        self.air.n=self.V_Cylinder/self.air.State.v  # calcualte number of moles of air
        self.air.m=self.air.n*self.air.MW

        self.State1=self.air.set(P=self.p_initial, T=self.T_initial)
        self.State2=self.air.set(v=self.State1.v/self.Ratio, s=self.State1.s)
        self.State3=self.air.set(T=self.T_high, v=self.State2.v)
        self.State4=self.air.set(v=self.State1.v, s=self.State3.s)
        
        self.W_Compression = self.air.n*(self.State2.u-self.State1.u)
        self.W_Power = self.air.n*(self.State3.u-self.State3.u)
        self.Q_In = self.air.n*(self.State3.u-self.State2.u)
        self.Q_Out = self.air.n*(self.State4.u-self.State1.u)
        
        self.W_Cycle=self.W_Power-self.W_Compression
        self.Eff=self.W_Cycle/self.Q_In

        self.upperCurve=StateDataForPlotting()
        self.lowerCurve=StateDataForPlotting()
        self.calculated=False
        self.cycleType='otto'

    def getSI(self):
        return self.units.SI
    
class ottoCycleController():
    def __init__(self, model=None, ax=None):
        self.model=ottoCycleModel() if model is None else model
        self.view=ottoCycleView()
        self.view.ax = ax

    #region Functions that operate on the model (i.e., change model state)
    def calc(self):
        # read values from the GUI
        T0=float(self.view.le_TLow.text())
        P0=float(self.view.le_P0.text())
        V0=float(self.view.le_V0.text())
        TH=float(self.view.le_THigh.text())
        CR=float(self.view.le_CR.text())
        metric=self.view.rdo_Metric.isChecked()
        self.set(T_0=T0, P_0=P0, V_0=V0, T_High=TH, ratio=CR, SI=metric)

    def set(self, T_0=25.0, P_0=100.0, V_0=1.0, T_High=1500.0, ratio=6.0, SI=True):
        """
        Sets the initial state of the air and converts units from input
        :param T_0: Initial temperature in absolute units (R or K)
        :param P_0: Initial pressure in (atm or pa)
        :param V_0: Initial volume in (ft^3 or m^3)
        :param T_High: High temperature in (R or K)
        :param ratio: Compression ratio
        :param SI: boolean
        :return: none
        """
        self.model.units.set(SI=SI)
        self.model.T_initial=T_0 if SI else T_0/self.model.units.CF_T
        self.model.p_initial=P_0 if SI else P_0/self.model.units.CF_P
        self.model.T_high=T_High if SI else T_High/self.model.units.CF_T
        self.model.V_Cylinder=V_0 if SI else V_0/self.model.units.CF_V
        self.model.Ratio=ratio

        #note that all state calculations are for molar values
        self.model.State1=self.model.air.set(P=self.model.p_initial, T=self.model.T_initial, name='State 1 - BDC')
        self.model.State2=self.model.air.set(v=self.model.State1.v/self.model.Ratio, s=self.model.State1.s, name='State 2 - TDC')
        self.model.State3=self.model.air.set(T=self.model.T_high, v=self.model.State2.v, name='State 3 - TDC')
        self.model.State4=self.model.air.set(v=self.model.State1.v, s=self.model.State3.s, name='State 4 - BDC')

        self.model.air.n=self.model.V_Cylinder/self.model.air.State.v  # calcualte number of moles of air
        self.model.air.m=self.model.air.n*self.model.air.MW

        self.model.W_Compression = self.model.State2.u - self.model.State1.u
        self.model.W_Power = self.model.State3.u - self.model.State4.u
        self.model.Q_In = self.model.State3.u - self.model.State2.u
        self.model.Q_Out = self.model.State4.u - self.model.State1.u

        self.model.W_Cycle = self.model.W_Power - self.model.W_Compression
        self.model.Eff = 100.0*self.model.W_Cycle / self.model.Q_In
        self.model.calculated=True

        self.buildDataForPlotting()
        self.updateView()

    def buildDataForPlotting(self ):
        """
        I want to create state data between states 1-2, 2-3, 3-4, 4-1
        I'll piece together an upperCurve data set from 2-3, 3-4, 4-1
        The lowerCurve data set is 1-2
        :return:
        """
        # clear out any old data
        self.model.upperCurve.clear()
        self.model.lowerCurve.clear()
        a = air()  # an air object
        #region build upperCurve
        # region states from 2-3 (v=const, T from T2->T3)
        DeltaT=np.linspace(self.model.State2.T, self.model.State3.T, 30)
        for T in DeltaT:
            state = a.set(T=T, v=self.model.State2.v)
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))
        # endregion
        # region states from 3-4 (v=from TDC to BDC, s=const.)
        DeltaV = np.linspace(self.model.State3.v, self.model.State4.v, 30)
        for v in DeltaV:
            state = a.set(v=v, s=self.model.State3.s)
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))
        # endregion
        # region states from 4-1 (v=const, T from T4->T1)
        DeltaT=np.linspace(self.model.State4.T, self.model.State1.T, 30)
        for T in DeltaT:
            state = a.set(T=T, v=self.model.State4.v)
            self.model.upperCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))
        # endregion
        #endregion

        #region build lowerCurve
        # region states from 1-2 (v=from BDC to TDC, s=const.)
        DeltaV=np.linspace(self.model.State1.v, self.model.State2.v, 30)
        for v in DeltaV:
            state = a.set(v=v, s=self.model.State1.s)
            self.model.lowerCurve.add((state.T, state.P, state.u, state.h, state.s, state.v))
        # endregion
        #endregion
    #endregion

    # region Functions that operate on the view
    def plot_cycle_XY(self, X='s', Y='T', logx=False, logy=False, mass=False, total=False):
        self.view.plot_cycle_XY(self.model, X=X, Y=Y, logx=logx,logy=logy, mass=mass, total=total)
    
    def print_summary(self):
        self.view.print_summary(self.model)
    
    def get_summary(self):
        return self.view.get_summary(self.model)

    def setWidgets(self, w=None):
        tlow=w[6].text()
        [self.view.lbl_THigh, self.view.lbl_TLow, self.view.lbl_P0, self.view.lbl_V0, self.view.lbl_CR,
        self.view.le_THigh, self.view.le_TLow, self.view.le_P0, self.view.le_V0, self.view.le_CR,
        self.view.le_T1, self.view.le_T2, self.view.le_T3, self.view.le_T4,
        self.view.lbl_T1Units, self.view.lbl_T2Units, self.view.lbl_T3Units, self.view.lbl_T4Units,
        self.view.le_PowerStroke, self.view.le_CompressionStroke, self.view.le_HeatAdded, self.view.le_Efficiency,
        self.view.lbl_PowerStrokeUnits, self.view.lbl_CompressionStrokeUnits, self.view.lbl_HeatInUnits,
        self.view.rdo_Metric, self.view.cmb_Abcissa, self.view.cmb_Ordinate,
        self.view.chk_LogAbcissa, self.view.chk_LogOrdinate, self.view.ax, self.view.canvas]=w

        tlow=self.view.le_TLow.text()
        pass

    def updateView(self):
        self.view.updateView(cycle=self.model)
    #endregion

class ottoCycleView():
    def __init__(self):
        #region define some widgets
        self.lbl_THigh = qtw.QLabel()
        self.lbl_TLow = qtw.QLabel() 
        self.lbl_P0 = qtw.QLabel() 
        self.lbl_V0 = qtw.QLabel() 
        self.lbl_CR = qtw.QLabel()
        self.le_THigh = qtw.QLineEdit() 
        self.le_TLow = qtw.QLineEdit() 
        self.le_P0 = qtw.QLineEdit() 
        self.le_V0 = qtw.QLineEdit() 
        self.le_CR = qtw.QLineEdit()
        self.le_T1 = qtw.QLineEdit() 
        self.le_T2 = qtw.QLineEdit() 
        self.le_T3 = qtw.QLineEdit() 
        self.le_T4 = qtw.QLineEdit()
        self.lbl_T1Units = qtw.QLabel() 
        self.lbl_T2Units = qtw.QLabel() 
        self.lbl_T3Units = qtw.QLabel() 
        self.lbl_T4Units = qtw.QLabel()
        self.le_Efficiency = qtw.QLineEdit() 
        self.le_PowerStroke = qtw.QLineEdit()
        self.le_CompressionStroke=qtw.QLineEdit()
        self.le_HeatAdded = qtw.QLineEdit()
        self.lbl_PowerStrokeUnits=qtw.QLabel()
        self.lbl_CompressionStrokeUnits=qtw.QLabel()
        self.lbl_HeatInUnits = qtw.QLabel()
        self.rdo_Metric = qtw.QRadioButton()
        self.cmb_Abcissa = qtw.QComboBox()
        self.cmb_Ordinate = qtw.QComboBox()
        self.chk_LogAbcissa = qtw.QCheckBox()
        self.chk_LogOrdinate = qtw.QCheckBox()
        self.canvas=None
        self.ax=None
        #endregion

    def updateView(self, cycle):
        cycle.units.set(SI=self.rdo_Metric.isChecked())
        logx=self.chk_LogAbcissa.isChecked()
        logy=self.chk_LogOrdinate.isChecked()
        xvar=self.cmb_Abcissa.currentText()
        yvar=self.cmb_Ordinate.currentText()
        if cycle.calculated:
            self.plot_cycle_XY(cycle, X=xvar, Y=yvar, logx=logx, logy=logy, mass=False, total=True)
        self.updateDisplayWidgets(Model=cycle)

    def print_summary(self, cycle):
        print('Cycle Summary for: ', cycle.name)
        print('\tEfficiency: {:0.3f}%'.format(cycle.Eff))
        print('\tPower Stroke: {:0.3f} kJ/kg'.format(cycle.W_Power))
        print('\tCompression Stroke: {:0.3f} kJ/kg'.format(cycle.W_Compression))
        print('\tHeat Added: {:0.3f} kJ/kg'.format(cycle.Q_In))
        cycle.State1.print()
        cycle.State2.print()
        cycle.State3.print()
        cycle.State4.print()

    def get_summary(self, cycle):
        '''
        This returns a formatted string to put on the plot of the otto cycle.
        :return:
        '''
        s = r'Summary:'
        s += '\n$\eta$: {:0.1f}% '.format(cycle.efficiency)
        s += '\n$\eta_{turbine}$: ' + '{:0.2f}'.format(cycle.eff_turbine) if cycle.eff_turbine < 1.0 else ''
        s += '\n$W_{turbine}$: ' + '{:0.1f} kJ/k'.format(cycle.turbine_work)
        s += '\n$W_{pump}$: ' + '{:0.1f} kJ/kg'.format(cycle.pump_work)
        s += '\n$Q_{boiler}$: ' + '{:0.1f} kJ/kg'.format(cycle.heat_added)
        return s

    def convertDataCol(self, cycle, data=None, colName='T', mass=False, total=False):
        UC=cycle.units
        n=cycle.air.n
        MW=cycle.air.MW
        TCF = 1.0 if UC.SI else UC.CF_T
        PCF = 1.0 if UC.SI else UC.CF_P
        hCF = 1.0 if UC.SI else UC.CF_e
        uCF = 1.0 if UC.SI else UC.CF_e
        sCF = 1.0 if UC.SI else UC.CF_s
        vCF = 1.0 if UC.SI else UC.CF_v  # convert m^3/mol to ft^3/lbmol
        nCF = 1.0 if UC.SI else UC.CF_n  # convert mol to lbmol
        mCF = 1.0 if UC.SI else UC.CF_Mass
        if mass:
            hCF/=MW  # kJ/kmol to kJ/kg  or Btu/lbmol to Btu/lbmass
            uCF/=MW
            sCF/=MW
            vCF/=MW
        elif total:
            hCF*=n*nCF
            uCF*=n*nCF
            sCF*=n*nCF
            vCF*=n*nCF
        w=colName.lower()
        if w=='t':
            return [T*TCF for T in data]
        if w=='h':
            return [h*hCF for h in data]
        if w=='u':
            return [u*uCF for u in data]
        if w=='s':
            return [s*sCF for s in data]
        if w=='v':
            return [v*vCF for v in data]
        if w=='p':
            return [P*PCF for P in data]

    def plot_cycle_XY(self, cycle, X='s', Y='T',logx=False, logy=False, mass=False, total=False):
        """
        I want to plot any two thermodynaimc properties on X and Y
        Data is in molar metric units.  I may need to convert it.
        :param X: letter for which variable to plot on X axis
        :param Y: letter for which variable to plot on Y axis
        :return:
        """
        if X==Y:
            return
        QTPlotting = True  # assumes we are plotting onto a QT GUI form
        if self.ax == None:
            self.ax = plt.subplot()
            QTPlotting = False  # actually, we are just using CLI and showing the plot

        ax = self.ax
        ax.clear()

        ax.set_xscale('log' if logx else 'linear')
        ax.set_yscale('log' if logy else 'linear')

        # plot the upper and lower curves
        XdataLC=self.convertDataCol(cycle, colName=X,data=cycle.lowerCurve.getDataCol(X), mass=mass, total=total)
        YdataLC=self.convertDataCol(cycle, colName=Y,data=cycle.lowerCurve.getDataCol(Y), mass=mass, total=total)
        XdataUC=self.convertDataCol(cycle, colName=X,data=cycle.upperCurve.getDataCol(X), mass=mass, total=total)
        YdataUC=self.convertDataCol(cycle, colName=Y,data=cycle.upperCurve.getDataCol(Y), mass=mass, total=total)
        ax.plot(XdataLC, YdataLC, color='k')
        ax.plot(XdataUC, YdataUC, color='g')

        # add axis labels
        cycle.units.setPlotUnits(SI=cycle.units.SI, mass=mass, total=total)
        ax.set_ylabel(cycle.lowerCurve.getAxisLabel(Y, Units=cycle.units), fontsize='large')
        ax.set_xlabel(cycle.lowerCurve.getAxisLabel(X, Units=cycle.units), fontsize='large')

        # put a title on the plot
        cycle.name = 'Otto Cycle'
        ax.set_title(cycle.name, fontsize='large')

        # modify the tick marks
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True,labelsize='large')

        # plot the circles for states 1, 2, 3, and 4
        state1=dc(cycle.State1)
        state1.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)
        state2=dc(cycle.State2)
        state2.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)
        state3=dc(cycle.State3)
        state3.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)
        state4=dc(cycle.State4)
        state4.ConvertStateData(SI=cycle.getSI(), Units=cycle.units, n=cycle.air.n, MW=cycle.air.MW, mass=mass, total=total)

        ax.plot(state1.getVal(X), state1.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(state2.getVal(X), state2.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(state3.getVal(X), state3.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')
        ax.plot(state4.getVal(X), state4.getVal(Y), marker='o', markerfacecolor='w', markeredgecolor='k')
        # # set limits on x and y
        xmin = min(min(cycle.upperCurve.getDataCol(X)), min(cycle.lowerCurve.getDataCol(X)))
        xmax = max(max(cycle.upperCurve.getDataCol(X)), max(cycle.lowerCurve.getDataCol(X)))
        ymin=min(min(cycle.upperCurve.getDataCol(Y)), min(cycle.lowerCurve.getDataCol(Y)))
        ymax=max(max(cycle.upperCurve.getDataCol(Y)), max(cycle.lowerCurve.getDataCol(Y)))
        #ax.set_xlim(xmin,xmax)
        #ax.set_ylim(ymin,ymax)
        deltax=xmax-xmin
        deltay=ymax-ymin
        # add the summary text to the plot
        # ax.text(xmin+0.05*deltax, ymin+0.7*deltay, self.get_summary(cycle))
        # show the plot
        if QTPlotting == False:
            plt.show()
        else:
            self.canvas.draw()

    def updateDisplayWidgets(self, Model=None):
        # fill out the temperature values
        U=Model.units
        SI=U.SI

        self.lbl_THigh.setText('T High ({})'.format(Model.units.TUnits))
        self.lbl_TLow.setText('T Low ({})'.format(Model.units.TUnits))
        self.lbl_P0.setText('P0 ({})'.format(Model.units.PUnits))
        self.lbl_V0.setText('V0 ({})'.format(Model.units.VUnits))

        self.lbl_T1Units.setText(Model.units.TUnits)
        self.lbl_T2Units.setText(Model.units.TUnits)
        self.lbl_T3Units.setText(Model.units.TUnits)
        self.lbl_T4Units.setText(Model.units.TUnits)

        if Model.units.changed or Model.calculated:
            if Model.calculated:
                CFE = 1.0 if SI else U.CF_E
                CFP = 1.0 if SI else U.CF_P
                CFV = 1.0 if SI else U.CF_V
                self.le_TLow.setText(('{:0.2f}'.format(Model.T_initial if SI else U.T_KtoR(Model.T_initial))))
                self.le_THigh.setText(('{:0.2f}'.format(Model.T_high if SI else U.T_KtoR(Model.T_high))))
                self.le_P0.setText('{:0.2f}'.format(Model.p_initial * CFP))
                self.le_V0.setText('{:0.4f}'.format(Model.V_Cylinder*CFV))

                self.le_T1.setText('{:0.2f}'.format(Model.State1.T if SI else U.T_KtoR(Model.State1.T)))
                self.le_T2.setText('{:0.2f}'.format(Model.State2.T if SI else U.T_KtoR(Model.State2.T)))
                self.le_T3.setText('{:0.2f}'.format(Model.State3.T if SI else U.T_KtoR(Model.State3.T)))
                self.le_T4.setText('{:0.2f}'.format(Model.State4.T if SI else U.T_KtoR(Model.State4.T)))

                # fill out the other properties for the diesel cycle
                self.le_Efficiency.setText('{:0.3f}'.format(Model.Eff))
                self.le_PowerStroke.setText('{:0.3f}'.format(Model.air.n*Model.W_Power*CFE))
                self.le_CompressionStroke.setText('{:0.3f}'.format(Model.air.n*Model.W_Compression*CFE))
                self.le_HeatAdded.setText('{:0.3f}'.format(Model.air.n*Model.Q_In*CFE))
                self.lbl_PowerStrokeUnits.setText(Model.units.EUnits)
                self.lbl_CompressionStrokeUnits.setText(Model.units.EUnits)
                self.lbl_HeatInUnits.setText(Model.units.EUnits)
            else:
                CFE = 1/U.CF_E if SI else U.CF_E
                CFP = 1/U.CF_P if SI else U.CF_P
                CFV = 1/U.CF_V if SI else U.CF_V
                t_high=float(self.le_THigh.text())
                t_initial=float(self.le_TLow.text())
                p_initial=float(self.le_P0.text())
                v_initial=float(self.le_V0.text())
                self.le_THigh.setText(('{:0.2f}'.format(U.T_RtoK(t_high) if SI else U.T_KtoR(t_high))))
                self.le_TLow.setText(('{:0.2f}'.format(U.T_RtoK(t_initial) if SI else U.T_KtoR(t_initial))))
                self.le_P0.setText('{:0.2f}'.format(p_initial * CFP))
                self.le_V0.setText('{:0.4f}'.format(v_initial * CFV))
            Model.units.changed=False
#endregion

def main():
    oc=ottoCycleController()
    oc.set(T_0=540.0, P_0=1.0, T_High=3600.0, ratio=8.0, V_0=0.02, SI=False)
    oc.plot_cycle_XY(X='v', Y='P', total=True)
    #oc.print_summary()
    #oc.plot_cycle_XY(X='s', Y='T', mass=True)

if __name__ == "__main__":
    app = qtw.QApplication(sys.argv)
    main()
