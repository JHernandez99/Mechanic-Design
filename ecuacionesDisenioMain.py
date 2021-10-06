from PyQt5 import QtWidgets, uic,QtGui
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt,QThread, pyqtSignal
import sys
import numpy as np
pi = 3.14159265359

class Ui(QtWidgets.QMainWindow):
    def __init__(self):#Constructor
        #Datos para cargar la interfaz y abrirla
        super(Ui,self).__init__()
        self.ecuacionesDisenioGUI = uic.loadUi("mainGUI.ui",self)
        self.ecuacionesDisenioGUI.show()
        self.setWindowTitle("Diseño Mecánico - Ecuaciones")
        #Instrucciones para llamar a los widgets
        #self.numero0 = self.findChild(QtWidgets.QPushButton, 'numero0')
        #self.numero0.clicked.connect(self.numero0Presionado)

        #display
        #self.displayLCD = self.findChild(QtWidgets.QLCDNumber,'displayLCD')

        #boton de salir
        self.calcular = self.findChild(QtWidgets.QPushButton, 'btnCalcular')
        self.calcular.clicked.connect(self.calcularEcuaciones)
        self.borrarTodo = self.findChild(QtWidgets.QPushButton, 'btnBorrarTodo')
        self.borrarTodo.clicked.connect(self.borrarTodoFun)


    def borrarTodoFun(self):

        pass

    def calcularEcuaciones(self):
        self.tblDatos = self.findChild(QtWidgets.QTableWidget, 'tblDatos')
        self.tblResultados = self.findChild(QtWidgets.QTableWidget, 'tblResultados')


        #print(self.tblDatos.item(1,1).text())
        SutPa = float(self.tblDatos.item(0,1).text())

        SutPsi = float(self.tblDatos.item(0,2).text())
        SytPa = float(self.tblDatos.item(1, 1).text())
        SytPsi = float(self.tblDatos.item(1, 2).text())
        SePa = float(self.tblDatos.item(2, 1).text())
        SePsi = float(self.tblDatos.item(2, 2).text())
        MmPa = float(self.tblDatos.item(4, 1).text())
        MmPsi = float(self.tblDatos.item(4, 2).text())
        MaPa = float(self.tblDatos.item(5, 1).text())
        MaPsi = float(self.tblDatos.item(5, 2).text())
        TmPa = float(self.tblDatos.item(6, 1).text())
        TmPsi = float(self.tblDatos.item(6, 2).text())
        TaPa = float(self.tblDatos.item(7, 1).text())
        TaPsi = float(self.tblDatos.item(7, 2).text())
        KfPa = float(self.tblDatos.item(9, 1).text())
        KfPsi = float(self.tblDatos.item(9, 2).text())
        KfsPa = float(self.tblDatos.item(10, 1).text())
        KfsPsi = float(self.tblDatos.item(10, 2).text())
        dM = float(self.tblDatos.item(12, 1).text())
        dIn = float(self.tblDatos.item(12, 2).text())
        n = float(self.tblDatos.item(14, 2).text())

        #ED- GOODMAN En Pascales
        EDGoodmanInv = (16/(pi*(dM**3)))*(((1/SePa)*(4*(KfPa*MaPa)**2 + 3*(KfsPa*MmPa)**2)**0.5) + ((1/SutPa)*(4*(KfPa*MmPa)**2 + 3*(KfsPa*TmPa)**2)**0.5))
        EDGoodmanPa = 1/EDGoodmanInv

        EDGoodmanInvPsi = (16 / (pi * (dIn ** 3))) * (
                    ((1 / SePsi) * (4 * (KfPsi * MaPsi) ** 2 + 3 * (KfsPsi * MmPsi) ** 2) ** 0.5) + (
                        (1 / SutPsi) * (4 * (KfPsi * MmPsi) ** 2 + 3 * (KfsPsi * TmPsi) ** 2) ** 0.5))
        EDGoodmanPsi = 1 / EDGoodmanInvPsi

        dGoodmanM = ((16*n)/(pi))*(((1/SePa)*(4*(KfPa*MaPa)**2 + 3*(KfsPa*MmPa)**2)**0.5) + ((1/SutPa)*(4*(KfPa*MmPa)**2 + 3*(KfsPa*TmPa)**2)**0.5))
        dGoodmanMi = dGoodmanM**(1/3)
        print(dGoodmanM**(1/3))
        dGoodmanIn = (16*n / (pi)) * (((1 / SePsi) * (4 * (KfPsi * MaPsi) ** 2 + 3 * (KfsPsi * MmPsi) ** 2) ** 0.5) + ((1 / SutPsi) * (4 * (KfPsi * MmPsi) ** 2 + 3 * (KfsPsi * TmPsi) ** 2) ** 0.5))
        dGoodmanIn = dGoodmanIn**(1/3)
        print(dGoodmanIn)

        #ED-Gerber
        A = np.sqrt(4*((KfPa*MaPa)**2) + 3*((KfsPa*TaPa)**2))
        B = np.sqrt(4*((KfPa*MmPa)**2) + 3*((KfsPa*TmPa)**2))
        EDGerberInv =((8*A)/(pi*(dM**3)*SePa))*( 1 + (( 1 + (((2*B*SePa)/(A*SutPa))**2))**(1/2)) )
        EDGerberPa = 1 / EDGerberInv

        dGerberM = ((8*n*A)/(pi*SePa))*( 1 + (( 1 + (((2*B*SePa)/(A*SutPa))**2))**(1/2)) )
        dGerberMi = dGerberM**(1/3)

        APsi = np.sqrt(4 * (KfPsi * MaPsi) ** 2 + 3 * (KfsPsi * TaPsi) ** 2)
        BPsi = np.sqrt(4 * (KfPsi * MmPsi) ** 2 + 3 * (KfsPsi * TmPsi) ** 2)
        EDGerberInvPsi = ((8 * APsi) / (pi * (dIn ** 3) * SePsi)) * (1 + (1 + (((2 * BPsi * SePsi) / (APsi * SutPsi)) ** 2)) ** 0.5)
        EDGerberPsi = 1 / EDGerberInvPsi

        dGerberIn =((8 * n* APsi) / (pi * SePsi)) * (1 + (1 + (((2 * BPsi * SePsi) / (APsi * SutPsi)) ** 2)) ** 0.5)
        dGerberIn = dGerberIn**(1/3)

        #ED-ASME eliptica
        EDASMEElipticaInv = (16/pi*(dM**3))*((4*((KfPa*MaPa)/SePa)**2 + 3*((KfsPa*TaPa)/SePa)**2 + 4*((KfPa*MmPa)/SytPa)**2 + 3*((KfsPa*TmPa)/SytPa)**2  )**0.5)
        EDASMEElipticaPa = 1/ EDASMEElipticaInv

        dASMEElipticaM = (16*n/pi)*((4*((KfPa*MaPa)/SePa)**2 + 3*((KfsPa*TaPa)/SePa)**2 + 4*((KfPa*MmPa)/SytPa)**2 + 3*((KfsPa*TmPa)/SytPa)**2  )**0.5)
        dASMEElipticaMi = dASMEElipticaM**(1/3)
        EDASMEElipticaInvPsi = (16 / pi * (dIn ** 3)) * ((4 * ((KfPsi * MaPsi) / SePsi) ** 2 + 3 * (
                    (KfsPsi * TaPsi) / SePsi) ** 2 + 4 * ((KfPsi * MmPsi) / SytPsi) ** 2 + 3 * (
                                                                  (KfsPsi * TmPsi) / SytPsi) ** 2) ** 0.5)
        EDASMEElipticaPsi = 1 / EDASMEElipticaInvPsi

        dASMEElipticaIn =(16*n / pi) * ((4 * ((KfPsi * MaPsi) / SePsi) ** 2 + 3 * (
                    (KfsPsi * TaPsi) / SePsi) ** 2 + 4 * ((KfPsi * MmPsi) / SytPsi) ** 2 + 3 * (
                                                                  (KfsPsi * TmPsi) / SytPsi) ** 2) ** 0.5)
        dASMEElipticaIn= dASMEElipticaIn**(1/3)

        #ED-Soderberg
        EDSoderbergInv = (16/pi*(dM**3))*(((1/SePa) * ((4*(KfPa*MaPa)**2 + 3*(KfsPa*TaPa)**2)**0.5)) + ((1/SytPa)*( 4*(KfPa*MmPa)**2 + 3*(KfsPa*TmPa)**2)**0.5))
        EDSoderbergPa = 1 / EDSoderbergInv

        dSoderbergM = (16*n/pi)*(((1/SePa) * ((4*(KfPa*MaPa)**2 + 3*(KfsPa*TaPa)**2)**0.5)) + ((1/SytPa)*( 4*(KfPa*MmPa)**2 + 3*(KfsPa*TmPa)**2)**0.5))
        dSoderbergMi = dSoderbergM**(1/3)

        EDSoderbergInvPsi = (16 / pi * (dIn ** 3)) * (
                    ((1 / SePsi) * ((4 * (KfPsi * MaPsi) ** 2 + 3 * (KfsPsi * TaPsi) ** 2) ** 0.5)) + (
                        (1 / SytPsi) * (4 * (KfPsi * MmPsi) ** 2 + 3 * (KfsPsi * TmPsi) ** 2) ** 0.5))
        EDSoderbergPsi = 1 / EDSoderbergInvPsi

        dSoderbergIn = (16*n / pi) * (
                    ((1 / SePsi) * ((4 * (KfPsi * MaPsi) ** 2 + 3 * (KfsPsi * TaPsi) ** 2) ** 0.5)) + (
                        (1 / SytPsi) * (4 * (KfPsi * MmPsi) ** 2 + 3 * (KfsPsi * TmPsi) ** 2) ** 0.5))
        dSoderbergIn = dSoderbergIn**(1/3)

        #Von Mises
        VonMisesInv = (((32*KfPa*(MmPa + MaPa))/(pi*(dM**3)))**2 + 3*((16*KfsPa*(TmPa + TaPa))/(pi*(dM**3)))**2)**0.5
        VonMisesPa = SytPa/VonMisesInv

        dVonMisesMi = (((32*KfPa*(MmPa + MaPa))/(pi*(VonMisesPa)))**2 + 3*((16*KfsPa*(TmPa + TaPa))/(pi*(VonMisesPa)))**2)**(1/6)

        VonMisesInvPsi = (((32 * KfPsi * (MmPsi + MaPsi)) / (pi * (dIn ** 3))) ** 2 + 3 * (
                    (16 * KfsPsi * (TmPsi + TaPsi)) / (pi * (dIn ** 3))) ** 2) ** 0.5
        VonMisesPsi = SytPsi / VonMisesInvPsi

        dVonMisesIn = (((32 * KfPsi * (MmPsi + MaPsi)) / (pi * (VonMisesPsi))) ** 2 + 3 * (
                    (16 * KfsPsi * (TmPsi + TaPsi)) / (pi * (VonMisesPsi))) ** 2) ** (1/6)

        self.tblResultados.item(0, 0).setText(str(EDGoodmanPa))
        self.tblResultados.item(1, 0).setText(str(EDGerberPa))
        self.tblResultados.item(2, 0).setText(str(EDASMEElipticaPa))
        self.tblResultados.item(3, 0).setText(str(EDSoderbergPa))
        self.tblResultados.item(4, 0).setText(str(VonMisesPa))

        self.tblResultados.item(0, 1).setText(str(EDGoodmanPsi))
        self.tblResultados.item(1, 1).setText(str(EDGerberPsi))
        self.tblResultados.item(2, 1).setText(str(EDASMEElipticaPsi))
        self.tblResultados.item(3, 1).setText(str(EDSoderbergPsi))
        self.tblResultados.item(4, 1).setText(str(VonMisesPsi))

        self.tblResultados.item(6, 0).setText(str(dGoodmanMi))
        self.tblResultados.item(7, 0).setText(str(dGerberMi))
        self.tblResultados.item(8, 0).setText(str(dASMEElipticaMi))
        self.tblResultados.item(9, 0).setText(str(dSoderbergMi))
        self.tblResultados.item(10, 0).setText(str(dVonMisesMi))

        self.tblResultados.item(6, 1).setText(str(dGoodmanIn))
        self.tblResultados.item(7, 1).setText(str(dGerberIn))
        self.tblResultados.item(8, 1).setText(str(dASMEElipticaIn))
        self.tblResultados.item(9, 1).setText(str(dSoderbergIn))
        self.tblResultados.item(10, 1).setText(str(dVonMisesIn))

        print("Operacion con exito")


app = QtWidgets.QApplication(sys.argv)
windows = Ui()
app.exec()
#wrong base class of toplevel widget main ui qdialog
#rrevisar en qt la clase principal de la interface ya que puede ser QMainWindow o QDialog