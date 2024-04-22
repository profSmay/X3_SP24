import math
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui



from copy import deepcopy as dc

#region special classes for circuit elements
# Copied this position class from pipe network problem
class Position():
    """
    I made this position for holding a position in 3D space (i.e., a point).  I've given it some ability to do
    vector arithmitic and vector algebra (i.e., a dot product).  I could have used a numpy array, but I wanted
    to create my own.  This class uses operator overloading as explained in the class.
    """
    def __init__(self, pos=None, x=None, y=None, z=None):
        """
        x, y, and z have the expected meanings
        :param pos: a tuple (x,y,z)
        :param x: float
        :param y: float
        :param z: float
        """
        # set default values
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        # unpack position from a tuple if given
        if pos is not None:
            self.x, self.y, self.z = pos
        # override the x,y,z defaults if they are given as arguments
        self.x = x if x is not None else self.x
        self.y = y if y is not None else self.y
        self.z = z if z is not None else self.z
        self.pt = qtc.QPointF(self.x, self.y)

    # region operator overloads $NEW$ 4/7/21
    # this is overloading the addition operator.  Allows me to add Position objects with simple math: c=a+b, where
    # a, b, and c are all position objects.
    def __add__(self, other):
        return Position((self.x + other.x, self.y + other.y, self.z + other.z))

    # this overloads the iterative add operator
    def __iadd__(self, other):
        if other in (float, int):
            self.x += other
            self.y += other
            self.z += other
            return self
        if type(other) == Position:
            self.x += other.x
            self.y += other.y
            self.z += other.z
            return self

    # this is overloading the subtract operator.  Allows me to subtract Positions. (i.e., c=b-a)
    def __sub__(self, other):
        return Position((self.x - other.x, self.y - other.y, self.z - other.z))

    # this overloads the iterative subtraction operator
    def __isub__(self, other):
        if other in (float, int):
            self.x -= other
            self.y -= other
            self.z -= other
            return self
        if type(other) == Position:
            self.x -= other.x
            self.y -= other.y
            self.z -= other.z
            return self

    # this is overloading the multiply operator.  Allows me to multiply a scalar or do a dot product (i.e., b=s*a or c=b*a)
    def __mul__(self, other):
        if type(other) in (float, int):
            return Position((self.x * other, self.y * other, self.z * other))
        if type(other) is Position:
            return Position((self.x * other.x, self.y * other.y, self.z * other.z))

    # this is overloading the __rmul__ operator so that s*Pt works.
    def __rmul__(self, other):
        return self * other

    # this is overloading the *= operator.  Same as a = Position((a.x*other, a.y*other, a.z*other))
    def __imul__(self, other):
        if type(other) in (float, int):
            self.x *= other
            self.y *= other
            self.z *= other
            return self

    # this is overloading the division operator.  Allows me to divide by a scalar (i.e., b=a/s)
    def __truediv__(self, other):
        if type(other) in (float, int):
            return Position((self.x / other, self.y / other, self.z / other))

    # this is overloading the /= operator.  Same as a = Position((a.x/other, a.y/other, a.z/other))
    def __idiv__(self, other):
        if type(other) in (float, int):
            self.x /= other
            self.y /= other
            self.z /= other
            return self

    def __round__(self, n=None):
        if n is not None:
            return Position(x=round(self.x, n), y=round(self.y, n), z=round(self.z, n))
        return self

    # endregion

    def set(self, strXYZ=None, tupXYZ=None, SI=True):
        """
            set position by string or tuple
        """
        lenCF = 1 if SI else 3.3
        if strXYZ is not None:
            cells = strXYZ.replace('(', '').replace(')', '').strip().split(',')
            x, y, z = float(cells[0]), float(cells[1]), float(cells[2])
            self.x = lenCF*float(x)
            self.y = lenCF*float(y)
            self.z = lenCF*float(z)
        elif tupXYZ is not None:
            x, y, z = tupXYZ  # [0], strXYZ[1],strXYZ[2]
            self.x = lenCF*float(x)
            self.y = lenCF*float(y)
            self.z = lenCF*float(z)

    def getTup(self):  # return (x,y,z) as a tuple
        return (self.x, self.y, self.z)

    def getStr(self, nPlaces=3, SI=True):
        lenCF=1 if SI else 3.3
        return '{}, {}, {}'.format(round(self.x*lenCF, nPlaces), round(self.y*lenCF, nPlaces), round(self.z*lenCF, nPlaces))

    def mag(self):  # normal way to calculate magnitude of a vector
        return (self.x ** 2 + self.y ** 2 + self.z ** 2) ** 0.5

    def normalize(self):  # typical way to normalize to a unit vector
        l = self.mag()
        if l <= 0.0:
            return
        self.__idiv__(l)

    def normalize2D(self):
        self.z = 0.0
        self.normalize()

    def getAngleRad(self):
        """
        Gets angle of position relative to an origin (0,0) in the x-y plane
        :return: angle in x-y plane in radians
        """
        l = self.mag()
        if l <= 0.0:
            return 0
        if self.y >= 0.0:
            return math.acos(self.x / l)
        return 2.0 * math.pi - math.acos(self.x / l)

    def getAngleDeg(self):
        """
        Gets angle of position relative to an origin (0,0) in the x-y plane
        :return: angle in x-y plane in degrees
        """
        return 180.0 / math.pi * self.getAngleRad()

    def midPt(self, p2=None):
        """
        find midpoint between self and p2.
        :param p2: a position
        :return: position in middle between self and p2
        """
        return Position(x=self.x+0.5*(p2.x-self.x), y=self.y+0.5*(p2.y-self.y), z=self.z+0.5*(p2.z-self.z))

class circuitNode():
    def __init__(self, name='none', x=0, y=0):
        self.name = name
        self.position = Position((x,y,0))
        self.draw = True

class resistor():
    def __init__(self, name='none', R=10, node1Name='a', node2Name='b'):
        self.name=name
        self.R=R
        self.node1Name=node1Name
        self.node2Name=node2Name

class inductor():
    def __init__(self, name='none', L=20, node1Name='a', node2Name='b'):
        self.name=name
        self.L=L
        self.node1Name=node1Name
        self.node2Name=node2Name

class QGraphicsArcItem(qtw.QGraphicsEllipseItem):
    def paint(self,painter, option,widget):
        """
        This overrides the paint function for a QGraphicsEllipseItem.
        The default paint for an ellipse actually draws a wedge rather than just the arc.
        :param painter:
        :param option:
        :param widget:
        :return:
        """
        painter.setPen(self.pen())
        painter.setBrush(self.brush())
        painter.drawArc(self.rect(), self.startAngle(), self.spanAngle())

class capacitor():
    def __init__(self, name='none', C=0.05, node1Name='a', node2Name='b'):
        self.name=name
        self.C=C
        self.node1Name=node1Name
        self.node2Name=node2Name

class wire():
    def __init__(self, name='none', node1Name='a', node2Name='b'):
        self.name=name
        self.node1Name = node1Name
        self.node2Name = node2Name

class voltageSource():
    def __init__(self, name='none', node1Name='a', node2Name='2'):
        self.name=name
        self.node1Name=node1Name
        self.node2Name=node2Name
#endreigon

#region circuit MVC
class circuitModel():
    def __init__(self):
        self.nodes=[]
        self.resistors=[]
        self.capacitors=[]
        self.inductors=[]
        self.voltageSources=[]
        self.wires=[]

    def getNode(self, name=None):
        for n in self.nodes:
            if n.name == name:
                return n

class circuitView():
    def __init__(self, dw=None):
        if dw is not None:
            self.setWidgets(dw)
        self.setupPensAndBrushes()
        self.setupScene()
        self.drawAGrid()
        self.scale=2.0
        self.Zoom=1.0
        self.elementFraction = 1/2
    def setupPensAndBrushes(self):
        # region setup pens and brushes and scene #copied from pipe network program
        # make the pens first
        # a thick darkGray pen
        self.penPipe = qtg.QPen(qtc.Qt.darkGray)
        self.penPipe.setWidth(4)
        # a fine black pen
        self.penNode = qtg.QPen(qtc.Qt.black)
        self.penNode.setStyle(qtc.Qt.SolidLine)
        self.penNode.setWidth(1)

        # a fine black pen
        self.penVS = qtg.QPen(qtc.Qt.black)
        self.penVS.setStyle(qtc.Qt.SolidLine)
        self.penVS.setWidth(1)

        # a medium black pen
        self.penConnect=qtg.QPen(qtc.Qt.black)
        self.penConnect.setStyle(qtc.Qt.SolidLine)
        self.penConnect.setWidth(2)

        # a pen for the grid lines
        self.penGridLines = qtg.QPen()
        self.penGridLines.setWidth(1)
        # I wanted to make the grid lines more subtle, so set alpha=25
        self.penGridLines.setColor(qtg.QColor.fromHsv(197, 144, 228, alpha=50))
        # now make some brushes
        # build a brush for filling with solid red
        self.brushFill = qtg.QBrush(qtc.Qt.darkRed)
        self.brushVS = qtg.QBrush(qtg.QColor.fromHsv(0,0,128, alpha=50))
        # a brush that makes a hatch pattern
        self.brushNode = qtg.QBrush(qtc.Qt.black)
        # a brush for the background of my grid
        self.brushGrid = qtg.QBrush(qtg.QColor.fromHsv(87, 98, 245, alpha=128))
        # Finally, the scene where pipe network objects are drawn
        self.scene = qtw.QGraphicsScene()
        # endregion

    def setWidgets(self, dw):
        self.gv_Main, self.lineEdit, self.form = dw

    def setupScene(self):
        # create a scene object
        self.scene = qtw.QGraphicsScene()
        self.scene.setObjectName("MyScene")
        self.scene.setSceneRect(-200, -200, 400, 400)  # xLeft, yTop, Width, Height
        self.scene.installEventFilter(self.form)  # install the event filter for use in mouse tracking
        self.gv_Main.setMouseTracking(True)

        # set the scene for the graphics view object
        self.gv_Main.setScene(self.scene)

    def setZoom(self):
        self.gv_Main.resetTransform()
        self.gv_Main.scale(self.Zoom, self.Zoom)

    def drawAGrid(self, DeltaX=5, DeltaY=5, Height=200, Width=200, CenterX=0, CenterY=0):
        """
        This makes a grid for reference.  No snapping to grid enabled.
        :param DeltaX: grid spacing in x direction
        :param DeltaY: grid spacing in y direction
        :param Height: height of grid (y)
        :param Width: width of grid (x)
        :param CenterX: center of grid (x, in scene coords)
        :param CenterY: center of grid (y, in scene coords)
        :param Pen: pen for grid lines
        :param Brush: brush for background
        :return: nothing
        """
        Pen = self.penGridLines
        Brush = self.brushGrid
        height = self.scene.sceneRect().height() if Height is None else Height
        width = self.scene.sceneRect().width() if Width is None else Width
        left = self.scene.sceneRect().left() if CenterX is None else (CenterX - width / 2.0)
        right = self.scene.sceneRect().right() if CenterX is None else (CenterX + width / 2.0)
        top = -1.0 * self.scene.sceneRect().top() if CenterY is None else (-CenterY + height / 2.0)
        bottom = -1.0 * self.scene.sceneRect().bottom() if CenterY is None else (-CenterY - height / 2.0)
        Dx = DeltaX
        Dy = DeltaY
        pen = qtg.QPen() if Pen is None else Pen

        # make the background rectangle first
        if Brush is not None:
            rect = qtw.QGraphicsRectItem(left, -top, width, height)
            rect.setBrush(Brush)
            rect.setPen(pen)
            self.scene.addItem(rect)
        # draw the vertical grid lines
        x = left
        while x <= right:
            lVert = qtw.QGraphicsLineItem(x, top, x, bottom)
            lVert.setPen(pen)
            self.scene.addItem(lVert)
            x += Dx
        # draw the horizontal grid lines
        y = bottom
        while y <= top:
            lHor = qtw.QGraphicsLineItem(left, -y, right, -y)  # now flip y
            lHor.setPen(pen)
            self.scene.addItem(lHor)
            y += Dy

    def drawCircuit(self, Model=circuitModel()):
        for n in Model.nodes:
            if n.draw:
                self.drawACircle(centerX=n.position.x*self.scale, centerY=n.position.y*self.scale, Radius=1*self.scale, brush=self.brushNode, pen=self.penNode, tooltip=n.name)
        for w in Model.wires:
            self.drawWire(w, Model=Model)
        for V in Model.voltageSources:
            self.drawVoltageSource(V,Model=Model)
        for R in Model.resistors:
            self.drawResistor(R,Model=Model)
        for C in Model.capacitors:
            self.drawCapacitor(C, Model=Model)
        for L in Model.inductors:
            self.drawInductor(L, Model=Model)

    def drawVoltageSource(self, V=voltageSource(), Model=circuitModel()):
        n1 = Model.getNode(V.node1Name)
        n2 = Model.getNode(V.node2Name)
        vec = n2.position-n1.position
        center = n1.position.midPt(n2.position)*self.scale
        dist=self.scale*vec.mag()
        radius=dist*self.elementFraction/2
        self.drawACircle(centerX=center.x , centerY=center.y , Radius=radius,
                         brush=self.brushVS, pen=self.penConnect, tooltip=V.name)
        #Draw +/- symbols
        self.drawALine(StartPt=Position(x=center.x-radius/3,y=center.y-radius/3), EndPt=Position(x=center.x+radius/3, y=center.y-radius/3), Pen=self.penConnect)
        self.drawALine(StartPt=Position(x=center.x-radius/3,y=center.y+radius/3), EndPt=Position(x=center.x+radius/3, y=center.y+radius/3), Pen=self.penConnect)
        self.drawALine(StartPt=Position(x=center.x,y=center.y), EndPt=Position(x=center.x, y=center.y+radius/3+radius/3), Pen=self.penConnect)

        self.drawALine(StartPt=n1.position*self.scale,EndPt=self.scale*(n1.position+(1-self.elementFraction)/2*vec), Pen=self.penConnect)
        self.drawALine(StartPt=n2.position*self.scale,EndPt=self.scale*(n2.position-(1-self.elementFraction)/2*vec), Pen=self.penConnect)

    def drawWire(self, w=wire(), Model=circuitModel()):
        n1 = Model.getNode(w.node1Name)
        n2 = Model.getNode(w.node2Name)
        self.drawALine(StartPt=n1.position*self.scale,EndPt=self.scale*(n2.position), Pen=self.penConnect)

    def drawResistor(self, R=resistor(), Model=circuitModel()):
        """
        For drawing a resistor symbol on a circuit diagram.
        :param R:
        :param Model:
        :return:
        """
        #$JES Missing Code$  draw a resistor
        pass

    def drawInductor(self, L=inductor(), Model=circuitModel()):
        """
        My code for drawing an inductor in a circuit diagram
        :param L:
        :param Model:
        :return:
        """
        #Step 1:  get node objects
        n1 = Model.getNode(L.node1Name)
        n2 = Model.getNode(L.node2Name)
        #Step 2:  calculate vector from node1 to node2
        vec=n2.position-n1.position
        #Step 3: determine orientation and containing rectangle of the inductor
        vert = True if n2.position.x == n1.position.x else False
        center = n1.position.midPt(n2.position)*self.scale
        dist=self.scale*vec.mag()
        width=self.elementFraction/2*dist if vert else self.elementFraction*dist
        height=self.elementFraction*dist if vert else self.elementFraction/2*dist
        #Step 4: draw the inductor symbol
        #region draw the inductor symbol
        if not vert:
            left = center.x - width / 2
            for i in range(4):
                if i == 0:
                    stAng=300
                    spnAng=180+60
                elif i == 3:
                    stAng=0
                    spnAng=180+60
                else:
                    stAng=300
                    spnAng=300
                self.drawAnArc(centerX=(left +(i+1)* width / 5), centerY=center.y, Radius=width/5, startAngle=stAng, spanAngle=spnAng, pen=self.penConnect)
        else:
            bottom=center.y-height/2
            for i in range(4):
                if i==0:
                    stAng=360-90
                    spnAng=180+60
                elif i == 3:
                    stAng=0-90-60
                    spnAng=180+60
                else:
                    stAng=300-90
                    spnAng=300
                self.drawAnArc(centerX=center.x, centerY=(bottom+(i+1)*(height/5)), Radius=height/5, startAngle=stAng, spanAngle=spnAng, pen=self.penConnect)
        #endregion
        #region connect inductor to nodes
        self.drawALine(StartPt=n1.position*self.scale,EndPt=self.scale*(n1.position+(1-self.elementFraction)/2*vec), Pen=self.penConnect)
        self.drawALine(StartPt=n2.position*self.scale,EndPt=self.scale*(n2.position-(1-self.elementFraction)/2*vec), Pen=self.penConnect)
        #endregion

    def drawCapacitor(self, C=capacitor(), Model=circuitModel()):
        #Step 1:  get node objects
        n1 = Model.getNode(C.node1Name)
        n2 = Model.getNode(C.node2Name)
        #Step 2:  calculate vector from node1 to node2
        vec=n2.position-n1.position
        #Step 3: determine orientation and containing rectangle of the inductor
        vert = True if n2.position.x == n1.position.x else False
        center = n1.position.midPt(n2.position)*self.scale
        dist=self.scale*vec.mag()
        width=self.elementFraction/2*dist if vert else self.elementFraction*dist
        height=self.elementFraction*dist if vert else self.elementFraction/2*dist
        #Step 4: draw the capacitor symbol
        #region draw the capacitor symbol
        gap=self.penConnect.width()*1
        if not vert:
            bottom = center.y - height / 2
            top = center.y+height/2
            self.drawALine(StartPt=Position(x=center.x-gap,y=bottom), EndPt=Position(x=center.x-gap, y=top), Pen=self.penConnect)
            self.drawALine(StartPt=Position(x=center.x+gap, y=bottom), EndPt=Position(x=center.x+gap,y=top), Pen=self.penConnect)
        else:
            left = center.x - width / 2
            right = center.x+width/2
            self.drawALine(StartPt=Position(x=left, y=center.y-gap), EndPt=Position(x=right, y=center.y-gap), Pen=self.penConnect)
            self.drawALine(StartPt=Position(x=left, y=center.y+gap), EndPt=Position(x=right, y=center.y+gap), Pen=self.penConnect)

        #endregion
        #region connect the capacitor symbol
        self.drawALine(StartPt=n1.position * self.scale,
                       EndPt=self.scale * (n1.position + (0.5-gap/dist) * vec), Pen=self.penConnect)
        self.drawALine(StartPt=n2.position * self.scale,
                       EndPt=self.scale * (n2.position - (0.5-gap/dist) * vec), Pen=self.penConnect)
        #endregion

    def drawACircle(self, centerX, centerY, Radius, startAngle=0, spanAngle=360, ccw=True, brush=None, pen=None, name=None, tooltip=None):
        scene = self.scene
        stAng = startAngle*16
        spnAngle = spanAngle*16
        spnAngle*=1 if ccw else -1

        # ellipse = qtw.QGraphicsEllipseItem(centerX - Radius, centerY - Radius, 2 * Radius, 2 * Radius)
        ellipse = qtw.QGraphicsEllipseItem(centerX - Radius, -1.0 * (centerY + Radius), 2 * Radius,
                                           2 * Radius)  # $NEW$ 4/7/21 flip y
        ellipse.setStartAngle(stAng)
        ellipse.setSpanAngle(spnAngle)

        if pen is not None:
            ellipse.setPen(pen)
        if brush is not None:
            ellipse.setBrush(brush)
        if name is not None:
            ellipse.setData(0, name)
        if tooltip is not None:
            ellipse.setToolTip(tooltip)
        scene.addItem(ellipse)

    def drawAnArc(self, centerX, centerY, Radius, startAngle=0, spanAngle=360, ccw=True, brush=None, pen=None, name=None, tooltip=None):
        scene = self.scene
        stAng = startAngle * 16
        spnAngle = spanAngle * 16
        spnAngle *= 1 if ccw else -1
        rect=qtc.QRectF(centerX - Radius, -1.0 * (centerY + Radius), 2 * Radius,2 * Radius)
        arc=QGraphicsArcItem()
        arc.setRect(rect)
        arc.setStartAngle(stAng)
        arc.setSpanAngle(spnAngle)

        if pen is not None:
            arc.setPen(pen)
        if brush is not None:
            arc.setBrush(brush)
        if name is not None:
            arc.setData(0, name)
        if tooltip is not None:
            arc.setToolTip(tooltip)
        scene.addItem(arc)

    def drawARectangle(self, centerX, centerY, height, width, brush=None, pen=None, name=None, tooltip=None):
        scene = self.scene
        # ellipse = qtw.QGraphicsEllipseItem(centerX - Radius, centerY - Radius, 2 * Radius, 2 * Radius)
        rect = qtw.QGraphicsRectItem(centerX - width/2, -1.0 * (centerY + height/2), width, height)  # flip y

        if pen is not None:
            rect.setPen(pen)
        if brush is not None:
            rect.setBrush(brush)
        if name is not None:
            rect.setData(0, name)
        if tooltip is not None:
            rect.setToolTip(tooltip)
        scene.addItem(rect)

    def drawALine(self, StartPt=Position(), EndPt=Position, Pen=qtg.QPen(qtg.QColor(qtc.Qt.black))):
        scene=self.scene
        ln = qtw.QGraphicsLineItem(StartPt.x, -StartPt.y, EndPt.x, -EndPt.y)
        ln.setPen(Pen)
        scene.addItem(ln)

class circuitController():
    def __init__(self, args):
        self.displayWidgets, self.inputWidgets = args
        self.spnd_Zoom = self.inputWidgets
        self.dlg =qtw.QFileDialog()

        self.Model = circuitModel()
        self.View = circuitView(self.displayWidgets)

    # region Functions responsible for importing/building objects from a data file
    def openFile(self):
        """
        Read the information from a circuit file.
        :return:
        """
        # open the file dialog box to search for the file I want
        filename = self.dlg.getOpenFileName(caption='Select a circuit file to open.')[0]
        if len(filename) == 0:  # no file selected
            return
        # self.le_FileName.setText(filename)  # echo the filename on the GUI
        file = open(filename, 'r')  # open the file
        data = file.readlines()  # read all the lines of the file into a list of strings
        self.importCircuit(data)  # import the pipe network information

    def importCircuit(self, data=[]):
        oldModel=dc(self.Model)
        self.Model=circuitModel()
        elementData=[]
        i=0
        while i < len(data):
            l=data[i].lower().strip()
            if l.find('<node')>=0:
                #region reading a node
                i+=1
                l=data[i]
                n=circuitNode()

                while l.lower().find('/node')<0:
                    if l.lower().find('name')>=0:
                        n.name=l.split(':')[1].strip()
                    elif l.lower().find('position')>=0:
                        cells=l.split(':')[1].split(',')
                        n.position.x=float(cells[0].strip())
                        n.position.y=float(cells[1].strip())
                    elif l.lower().find('draw')>=0:
                        TF= l.split(':')[1].strip().lower()
                        n.draw = TF.find('true') >= 0
                    i+=1
                    l=data[i]
                self.Model.nodes.append(dc(n))
                #endregion
            elif l.find('<resistor')>=0:
                #region reading a resistor
                #$JES MISSING CODE
                #endregion
            elif l.find('<inductor') >= 0:
                # region reading a inductor
                i += 1
                l = data[i]
                L = inductor()
                while l.lower().find('/inductor') < 0:
                    if l.lower().find('name') >= 0:
                        L.name = l.split(':')[1].strip()
                    elif l.lower().find('node1') >= 0:
                        L.node1Name = l.split(':')[1].strip()
                    elif l.lower().find('node2') >= 0:
                        L.node2Name = l.split(':')[1].strip()
                    i += 1
                    l = data[i]
                self.Model.inductors.append(dc(L))
                # endregion
            elif l.find('<capacitor') >= 0:
                # region reading a capacitor
                i += 1
                l = data[i]
                C = capacitor()
                while l.lower().find('/capacitor') < 0:
                    if l.lower().find('name') >= 0:
                        C.name = l.split(':')[1].strip()
                    elif l.lower().find('node1') >= 0:
                        C.node1Name = l.split(':')[1].strip()
                    elif l.lower().find('node2') >= 0:
                        C.node2Name = l.split(':')[1].strip()
                    i += 1
                    l = data[i]
                self.Model.capacitors.append(dc(C))
                # endregion
            elif l.find('<voltage source') >= 0:
                # region reading a voltage source
                i += 1
                l = data[i]
                V = voltageSource()
                while l.lower().find('/voltage source') < 0:
                    if l.lower().find('name') >= 0:
                        V.name = l.split(':')[1].strip()
                    elif l.lower().find('node1') >= 0:
                        V.node1Name = l.split(':')[1].strip()
                    elif l.lower().find('node2') >= 0:
                        V.node2Name = l.split(':')[1].strip()
                    i += 1
                    l = data[i]
                self.Model.voltageSources.append(dc(V))
                # endregion
            elif l.find('<wire') >= 0:
                # region reading a wire
                i += 1
                l = data[i]
                w = wire()
                while l.lower().find('/wire') < 0:
                    if l.lower().find('name') >= 0:
                        w.name = l.split(':')[1].strip()
                    elif l.lower().find('node1') >= 0:
                        w.node1Name = l.split(':')[1].strip()
                    elif l.lower().find('node2') >= 0:
                        w.node2Name = l.split(':')[1].strip()
                    i += 1
                    l = data[i]
                self.Model.wires.append(dc(w))
                # endregion
            i+=1
        self.drawCircuit()

    def setZoom(self):
        self.View.Zoom=self.spnd_Zoom.value()
        self.View.setZoom()

    def drawCircuit(self):
        self.View.drawCircuit(Model=self.Model)
    # endregion

    def getScene(self):
        return self.View.scene
#endregion
