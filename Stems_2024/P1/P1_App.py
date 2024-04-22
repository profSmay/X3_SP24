#region imports
from PyQt5.QtWidgets import QApplication, QWidget
import PyQt5.QtWidgets as qtw
import PyQt5.QtGui as qtg
from P1_GUI import Ui_MainForm  # you need to use pyuic5 to convert the .ui file to a .py file
from X2Q2_SP24 import simulate, doPlot
import sys
from Circuit_Classes import circuitController

#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
#endregion

#region class definitions
class main_window(Ui_MainForm, QWidget):
    def __init__(self):
        """
        Constructor for circuit simulator.
        """
        super().__init__()
        self.setupUi(self)
        # you should modify the window title appropriately
        #self.setWindowTitle("Circuit Simulator by Jim Smay (21 April, 2022)")

        self.inputWidgets = (self.le_Inductance, self.le_Resistance, self.le_Capacitence, self.le_Amplitude, self.le_Freq, self.le_Phase, self.le_simTime, self.le_simPts)
        self.displayWidgets = (self.layout_VertMain, self.layout_VertInput, self)
        self.controller = circuitController((self.inputWidgets,self.displayWidgets))
        self.setupSignalsAndSlots()
        self.show()


    def setupSignalsAndSlots(self):
        """
        Connect the push button to the calculate function.
        :return:
        """
        self.pb_Calculate.clicked.connect(self.calculate)
        pass

    def calculate(self):
        self.controller.calculate()

#endregion

if __name__ == "__main__":
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())