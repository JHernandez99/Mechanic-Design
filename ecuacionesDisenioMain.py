from PyQt5 import QtWidgets, uic
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
pi = 3.14159265359

class Ui(QtWidgets.QMainWindow):
    def __init__(self):#Constructor
        #Datos para cargar la interfaz y abrirla
        #self.fig, self.ax = plt.subplots(facecolor='gray')
        super(Ui,self).__init__()

        self.ecuacionesDisenioGUI = uic.loadUi("mainGUI.ui",self)
        self.ecuacionesDisenioGUI.show()
        self.setWindowTitle("Diseño Mecánico - Ecuaciones")
        self.calcular = self.findChild(QtWidgets.QPushButton, 'btnCalcular')
        self.calcular.clicked.connect(self.calcularEcuaciones)
        self.borrarTodo = self.findChild(QtWidgets.QPushButton, 'btnBorrarTodo')
        self.esfuerzoCortanteAdhesivo = self.findChild(QtWidgets.QPushButton, 'btnCalcular_2')
        self.esfuerzoCortanteAdhesivo.clicked.connect(self.calcularEsfuerzoCortanteAdhesivo)
        self.graficoInteractivo = self.findChild(QtWidgets.QPushButton, 'graficoInteractivo')
        self.graficoInteractivo.setEnabled(False)
        self.graficoInteractivo.clicked.connect(self.graficoInteractivoEsfuerzoCortanteAdhesivo)
        self.unionSoldadaFiletesTransversales = self.findChild(QtWidgets.QPushButton, 'btnCalcular_3')
        self.unionSoldadaFiletesTransversales.clicked.connect(self.calcularUnionSoldadaFiletesTransversales)
        self.graficoInteractivo3 = self.findChild(QtWidgets.QPushButton, 'graficoInteractivo_2')
        self.graficoInteractivo3.setEnabled(False)
        self.graficoInteractivo3.clicked.connect(self.graficoInteractivoUnionSoldadaFiletesTransversales)
    def calcularEcuaciones(self):
        self.tblDatos = self.findChild(QtWidgets.QTableWidget, 'tblDatos')
        self.tblResultados = self.findChild(QtWidgets.QTableWidget, 'tblResultados')

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
    def calcularEsfuerzoCortanteAdhesivo(self):

        self.tblDatos2 = self.findChild(QtWidgets.QTableWidget, 'tblDatos_2')
        G = float(self.tblDatos2.item(1,0).text())
        alfa = float(self.tblDatos2.item(2,1).text())
        h = float(self.tblDatos2.item(3,2).text())
        Eo = float(self.tblDatos2.item(5, 0).text())
        alfao = float(self.tblDatos2.item(6,1).text())
        to = float(self.tblDatos2.item(7, 2).text())
        Ei = float(self.tblDatos2.item(9,0).text())
        alfai = float(self.tblDatos2.item(10, 1).text())
        ti = float(self.tblDatos2.item(11, 2).text())
        l = float(self.tblDatos2.item(13,2).text())
        b = float(self.tblDatos2.item(14, 2).text())
        P = float(self.tblDatos2.item(15,0).text())
        deltaT = float(self.tblDatos2.item(17,0).text())
        w = np.sqrt((G/h)*( (1/(Eo*to)) + (2/(Ei*ti))))
        self.tblDatos2.item(16, 0).setText(str(w))
        self.x = np.array(np.arange(-0.5,0.5+0.001, 0.001))
        self.tth = []
        self.tp = []
        self.tcombinado = []
        contador = 0
        for i in self.x:
            self.tth.append(((alfai-alfao)*deltaT*w*np.sinh(w*i))/ ((1/(Eo*to) + 2/(Ei*ti))*np.cosh((w*l)/2)))
            self.tp.append((P*w*np.cosh(w*i))/(4*b*np.sinh((w*l)/2)))
            self.tcombinado.append(self.tth[contador]+self.tp[contador])
            contador+=1

        self.grafica = Canvas_grafica(self.x, self.tth,self.tp,self.tcombinado)

        gf = self.ecuacionesDisenioGUI.insertarGrafica.count()
        if gf >0:
            self.ecuacionesDisenioGUI.insertarGrafica.itemAt(0).widget().deleteLater()
        self.ecuacionesDisenioGUI.insertarGrafica.addWidget(self.grafica)
        self.graficoInteractivo.setEnabled(True)
    def graficoInteractivoEsfuerzoCortanteAdhesivo(self):

        plt.close("all")
        plt.plot(self.x,self.tth)
        plt.plot(self.x, self.tp)
        plt.plot(self.x, self.tcombinado)
        plt.grid()
        plt.show()
    def calcularUnionSoldadaFiletesTransversales(self):
        self.tblDatos3 = self.findChild(QtWidgets.QTableWidget, 'tblDatos_3')
        F = float(self.tblDatos3.item(1, 0).text())
        h = float(self.tblDatos3.item(2, 0).text())
        l = float(self.tblDatos3.item(3, 0).text())
        gradoInf = float(self.tblDatos3.item(5, 0).text())
        gradoSup = float(self.tblDatos3.item(6, 0).text())
        resolucion = float(self.tblDatos3.item(7, 0).text())
        self.theta  = np.array(np.arange(gradoInf, gradoSup+resolucion,resolucion))

        self.cortanteFileteTransversal = []
        self.normalFileteTransversal = []
        self.vonMisesFileteTransversal = []
        cuenta = 0
        for grado in self.theta:
            grado = grado*(np.pi/180)
            self.cortanteFileteTransversal.append((F/(h*l))*(np.sin(grado)*np.cos(grado) + (np.sin(grado)**2) ))
            self.normalFileteTransversal.append((F/(h*l))*((np.cos(grado)**2) +np.sin(grado) * np.cos(grado)))
            #self.vonMisesFileteTransversal.append(np.sqrt((F/(h*l))*(((np.cos(grado)**2) +np.sin(grado) * np.cos(grado))**2 + 3*(np.sin(grado)*np.cos(grado) + (np.sin(grado)**2) )**2 )))
            self.vonMisesFileteTransversal.append((self.normalFileteTransversal[cuenta]**2 + 3*self.cortanteFileteTransversal[cuenta]**2)**0.5)
            cuenta+=1
        grafica = Canvas_grafica_SoldadaFiletesTransversales(self.theta, self.cortanteFileteTransversal, self.normalFileteTransversal, self.vonMisesFileteTransversal )

        gf = self.ecuacionesDisenioGUI.insertarGrafica_2.count()
        if gf > 0:
            self.ecuacionesDisenioGUI.insertarGrafica_2.itemAt(0).widget().deleteLater()
        self.ecuacionesDisenioGUI.insertarGrafica_2.addWidget(grafica)
        self.graficoInteractivo3.setEnabled(True)
    def graficoInteractivoUnionSoldadaFiletesTransversales(self):
        print("abriendo")
        plt.close("all")
        plt.plot(self.theta, self.cortanteFileteTransversal)
        plt.plot(self.theta, self.normalFileteTransversal)
        plt.plot(self.theta, self.vonMisesFileteTransversal)
        plt.grid()
        plt.show()

class Canvas_grafica(FigureCanvasQTAgg):
    def __init__(self,x,tth,tp, tcombinado):

        self.fig, self.ax = plt.subplots(1,1)

        super().__init__(self.fig)
        self.ax.plot(x,tth)
        self.ax.plot(x, tp)
        self.ax.plot(x,tcombinado)
        self.ax.grid()
        self.ax.legend(["τth 'Termico'","τp 'Carga'","τc 'Combinado'"])

        self.fig.suptitle("Grafica del Esfuerzo Cortante en Union Traslapada con Adhesivo")
        self.ax.set(xlabel='X', ylabel='τ')
class Canvas_grafica_SoldadaFiletesTransversales(FigureCanvasQTAgg):
    def __init__(self,grados,cortanteFileteTransversal,normalFileteTransversal,vonMisesFileteTransversal ):
        self.fig, self.ax = plt.subplots(1,1)
        super().__init__(self.fig)
        self.ax.plot(grados, cortanteFileteTransversal)
        self.ax.plot(grados, normalFileteTransversal)
        self.ax.plot(grados, vonMisesFileteTransversal)
        self.ax.grid()
        self.ax.legend(['τ', 'σ', "σ'"])
        self.fig.suptitle("Grafica de Esfuerzos cortante, normal y Von Mises de una union soldada de filetes transversales")
        self.ax.set(xlabel='θ rad', ylabel="τ,σ,σ'")




app = QtWidgets.QApplication(sys.argv)
windows = Ui()
sys.exit(app.exec_())