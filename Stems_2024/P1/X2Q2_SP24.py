# region imports
from math import sin
import math
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
# endregion

# region function definitions
def odeSystem(t,X, *args):
    """
    this is the odeSystem callback I'm using for odeint().  The *args contains a callback for the input voltage.
    :param X: the current values of the state variables
    :param t: the current time from odeint
    :param args: fn(t), L, R, C
    :return: list of derivatives of state variables
    """
    #unpack X into convenient variables
    fn, L,R,C=args
    #assign friendly variable names
    i1=X[0]
    i2=X[1]
    #calculate the current input voltage
    vt=fn(t)
    #calculate derivatives for the state variables
    i1dot=#$JES MISSING CODE
    i2dot=#$JES MISSING CODE
    return [i1dot, i2dot]

def simulate(L=20, R=20,C=0.05, A=20, f=20, p=0, t=10, pts=500):
    """
    For simulating transient behavior of circuit.
    :param L: Inductance (H)
    :param R: Resistance (ohm)
    :param C: Capacitance (F)
    :param A: Amplitude (V)
    :param f: frequency (Hz)
    :param p: phase (deg)
    :param ax: axes for plotting
    :return: nothing
    """
    w = f*2*math.pi #frequency in radians/sec
    phi = p*math.pi/180.0 #phase in radians
    vin = lambda t: A*sin(w*t+phi)
    myargs=(vin, L, R, C)
    x0=[0,0] # initial values state variables
    tList=np.linspace(0,t,int(pts))  #time vector
    I=solve_ivp(odeSystem, t_span=[0,t], y0=x0, t_eval=tList, args=myargs)
    return I  #solve using solve_ivp

def doPlot(*args, ax=None):
    """
    Re-written on 4/21/2022 to adapt to plotting on GUI if ax is not None
    :param args: contains ((R, list of time values, and results of solve_ivp))
    :param ax:
    :return:
    """
    if ax == None:
        ax = plt.subplot()
        QTPlotting = False  # actually, we are just using CLI and showing the plot
    else:
        QTPlotting = True

    R,tList, I = args[0] # unpack results of simulation
    ax.clear()
    ax.plot(tList,I.y[0], linestyle='solid', color='k', label=r'$i_1(t)$')
    ax.plot(tList,I.y[1], linestyle='dashed', color='k', label=r'$i_2(t)$')
    ax.set_xlim(0,max(tList))
    minI=min(min(I.y[0]),min(I.y[1]))
    maxI=max(max(I.y[0]),max(I.y[1]))
    rangeI=abs(maxI-minI)
    ax.set_ylim(minI-0.01*rangeI, maxI+0.01*rangeI)
    ax.tick_params(axis='both', which='both', direction='in', top=True, labelsize=12)  # format tick marks
    ax.tick_params(axis='both', grid_linewidth=1, grid_linestyle='solid', grid_alpha=0.5)
    ax.tick_params(axis='both', which='minor')
    #ax.set_xticks([x for x in range(11)])
    ax.grid(True)
    ax.set_xlabel('t (s)', fontsize=12)
    ax.set_ylabel(r'$i_1, i_2 (A)$', fontsize=12)

    ax1=ax.twinx()
    yvals=R*(I.y[1]-I.y[0])
    yrange=abs(max(yvals)-min(yvals))
    ax1.plot(tList,yvals, linestyle='dotted', color='k', label=r'$v_c(t)$')
    ax1.set_ylim(min(yvals)-yrange*0.01,max(yvals)+yrange*0.01)
    ax1.tick_params(axis='y', which='both', direction='in', top=True, right=True, labelsize=12)  # format tick marks
    ax.legend(fontsize=12)
    ax1.legend(loc='lower right', fontsize=12)
    ax1.set_ylabel(r'$V_c(t) (V)$', fontsize=12)
    if not QTPlotting:
        plt.show()

def main():
    """
    For solving problem 2 on exam.  Note:  I'm passing a callback through solve_ivp to model in input voltage
    :return:
    """
    vin = lambda t:  20*sin(20*2.0*math.pi*t)  # a callback to be passed as an argument to odeSystem
    L=20
    R=10
    C=0.05
    myargs=(vin, L, R, C)  # a tuple containing:  (callback for voltage source, L, R, C)
    I0 = [0,0]  # initial conditions I[0]=0, I[1]=0
    tList=np.linspace(0,10,500)  # time vector
    I = simulate(L=L,R=R,C=C,A=20,f=20/(2*math.pi),p=0, t=10, pts=500)
    #I=solve_ivp(odeSystem,[0,10],I0,t_eval=tList, args=myargs)  #solve using solve_IVP
    VC = np.array(R*(I.y[1]-I.y[0])) # create a numpy array of voltage across capacitor as function of time
    #the plotting part
    doPlot((R,I.t, I))
    pass
# endregion

# region function calls
if __name__ ==  "__main__":
    main()
# endregion