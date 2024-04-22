from PyQt5.QtWidgets import QApplication, QWidget
import PyQt5.QtWidgets as qtw
import PyQt5.QtGui as qtg
import PyQt5.QtCore as qtc

from P2_GUI import Ui_Form
from Circuit_Classes import circuitController
import sys

class main_window(Ui_Form, QWidget):
    def __init__(self):
        """
        Constructor for circuit simulator.
        """
        #call constructor of parent class
        super().__init__()
        #setup UI from parent
        self.setupUi(self)

        self.displayWidgets = (self.gv_Main, self.lineEdit, self)
        self.inputWidgets = (self.spnd_Zoom)

        self.Controller=circuitController((self.displayWidgets, self.inputWidgets))
        self.setWindowTitle("Circuit Sketch by Jim Smay (21 April, 2022)")
        self.spnd_Zoom.valueChanged.connect(self.Controller.setZoom)
        self.pb_Open.clicked.connect(self.Controller.openFile)

        self.show()

    def eventFilter(self, obj, event):
        """
        This overrides the default eventFilter of the widget.  It takes action on events and then passes the event
        along to the parent widget.
        :param obj: The object on which the event happened
        :param event: The event itself
        :return: boolean from the parent widget
        """
        #region $NEW$ 4/6/21 for mouse tracking on the drawing
        if obj == self.Controller.View.scene:
            et=event.type()
            if event.type() == qtc.QEvent.GraphicsSceneMouseMove:
                w = app.topLevelAt(event.screenPos())
                scenePos=event.scenePos()
                s  =self.Controller.View.scene.itemAt(scenePos,self.gv_Main.transform())  # gets item from graphics scene under the mouse
                strScene="Mouse Position:  x = {}, y = {}".format(round(scenePos.x(),2), round(-scenePos.y(),2)) #$NEW$ 4/7/21 flip y
                if s is not None and s.data(0) is not None:  # when creating nodes and pipes, I used the setData() function to store a name
                    strScene += ' ' + s.data(0)
                self.lbl_MousePosition.setText(strScene)  # display information in a label
            if event.type() == qtc.QEvent.GraphicsSceneWheel:  # I added this to zoom on mouse wheel scroll
                if event.delta()>0:
                    self.spnd_Zoom.stepUp()
                else:
                    self.spnd_Zoom.stepDown()
                pass
        #endregion
        # allow table_Pipes, table_Nodes, tree_LoopPipes, and tree_Loops to respond to delete key
        if event.type() == qtc.QEvent.KeyPress:
            if event.key() == qtc.Qt.Key_Delete:
                if obj == self.table_Pipes:
                    self.deletePipe()
                elif obj == self.table_Nodes:
                    self.deleteNode()
                elif obj == self.tree_LoopPipes:
                    self.deleteLoopPipe()
                elif obj == self.tree_Loops:
                    self.deleteLoop()

        # pass the event along to the parent widget if there is one.
        return super(main_window, self).eventFilter(obj, event)


if __name__ == "__main__":
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())