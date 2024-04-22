# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'P2_GUI.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1260, 1183)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pb_Open = QtWidgets.QPushButton(Form)
        self.pb_Open.setObjectName("pb_Open")
        self.horizontalLayout.addWidget(self.pb_Open)
        self.lineEdit = QtWidgets.QLineEdit(Form)
        self.lineEdit.setObjectName("lineEdit")
        self.horizontalLayout.addWidget(self.lineEdit)
        self.spnd_Zoom = QtWidgets.QDoubleSpinBox(Form)
        self.spnd_Zoom.setSingleStep(0.1)
        self.spnd_Zoom.setProperty("value", 1.0)
        self.spnd_Zoom.setObjectName("spnd_Zoom")
        self.horizontalLayout.addWidget(self.spnd_Zoom)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.lbl_MousePosition = QtWidgets.QLabel(Form)
        self.lbl_MousePosition.setObjectName("lbl_MousePosition")
        self.verticalLayout.addWidget(self.lbl_MousePosition)
        self.gv_Main = QtWidgets.QGraphicsView(Form)
        self.gv_Main.setObjectName("gv_Main")
        self.verticalLayout.addWidget(self.gv_Main)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.pb_Open.setText(_translate("Form", "Open Circuit"))
        self.lbl_MousePosition.setText(_translate("Form", "Mouse Position"))

