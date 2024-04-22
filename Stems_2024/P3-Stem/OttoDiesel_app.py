#region imports
from OttoDiesel_GUI import Ui_Form
from PyQt5 import uic
import sys
from PyQt5 import QtWidgets as qtw
from Otto import ottoCycleController
from Diesel import dieselCycleController
from Air import *

#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
#endregion

class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """MainWindow constructor"""
        super().__init__()
        self.setupUi(self)
        # Main UI code goes here
        self.calculated=False

        #creating a canvas to draw a figure for the otto cycle
        self.figure=Figure(figsize=(8,8),tight_layout=True, frameon=True, facecolor='none')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.main_VerticalLayout.addWidget(self.canvas)

        #setting up some signals and slots
        self.rdo_Metric.toggled.connect(self.setUnits) #triggered when the state of the radio button changes
        self.btn_Calculate.clicked.connect(self.calcCycle)
        self.cmb_Abcissa.currentIndexChanged.connect(self.doPlot)
        self.cmb_Ordinate.currentIndexChanged.connect(self.doPlot)
        self.chk_LogAbcissa.stateChanged.connect(self.doPlot)
        self.chk_LogOrdinate.stateChanged.connect(self.doPlot)
        self.cmb_OttoDiesel.currentIndexChanged.connect(self.selectCycle)
        # End main ui code

        #create otto and diesel controller objects to work with later
        self.otto = #$JES MISSING CODE  # instantiate an ottoCycleController object
        self.diesel = #$JES MISSING CODE # instantiate a dieselCycleController object
        self.controller=self.otto
        self.someWidgets=[]

        self.someWidgets+=[self.lbl_THigh, self.lbl_TLow, self.lbl_P0, self.lbl_V0, self.lbl_CR]
        self.someWidgets+=[self.le_THigh, self.le_TLow, self.le_P0, self.le_V0, self.le_CR]
        self.someWidgets+=[self.le_T1, self.le_T2, self.le_T3, self.le_T4]
        self.someWidgets+=[self.lbl_T1Units, self.lbl_T2Units, self.lbl_T3Units, self.lbl_T4Units]
        self.someWidgets+=[self.le_PowerStroke, self.le_CompressionStroke, self.le_HeatAdded, self.le_Efficiency]
        self.someWidgets+=[self.lbl_PowerStrokeUnits, self.lbl_CompressionStrokeUnits, self.lbl_HeatInUnits]
        self.someWidgets+=[self.rdo_Metric, self.cmb_Abcissa, self.cmb_Ordinate]
        self.someWidgets+=[self.chk_LogAbcissa, self.chk_LogOrdinate, self.ax, self.canvas]
        #pass some widgets to the controller for both input and output
        self.otto.setWidgets(w=self.someWidgets)
        self.diesel.setWidgets(w=self.someWidgets)

        #show the form
        self.show()

    def clamp(self, val, low, high):
        if self.isfloat(val):
            val=float(val)
            if val>high:
                return float(high)
            if val <low:
                return float(low)
            return val
        return float(low)

    def isfloat(self,value):
        '''
        This function is a check to verify that a string can be converted to a float
        :return:
        '''
        if value=='NaN':return False
        try:
            float(value)
            return True
        except ValueError:
            return False

    def doPlot(self):
        self.controller.updateView()

    def selectCycle(self):
        otto = #$JES MISSING CODE # determine if otto cycle is chosen (true) or not (false -> diesel cycle)
        self.gb_Input.setTitle('Input for Air Standard {} Cycle:'.format('Otto' if otto else 'Diesel'))
        self.controller= #$JES MISSING CODE  # set self.controller to self.otto or self.diesel
        self.controller.updateView()

    def setUnits(self):
        self.controller.updateView()

    def calcCycle(self):
        '''
        This is called when the calculate button is clicked
        :return: nothing
        '''
        #calculate the cycle efficiency (and states 1,2,3,4)
        #$JES MISSING CODE call the calc function of the controller
        pass

#if this module is being imported, this won't run. If it is the main module, it will run.
if __name__== '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Otto Cycle Calculator')
    sys.exit(app.exec())