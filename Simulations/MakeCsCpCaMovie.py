import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod as TomasCircuit
from alive_progress import alive_bar

MHz_ = 1e6
pF_ = 1e-12

CsVals = np.linspace(25e-12,1000e-12,50)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,50)#capacitor is defined in steps
CaVals = np.linspace(100e-12,1000e-12,10)#capacitor is defined in steps

#leftmost = (np.abs(CpVals-220e-12)).argmin()
#rightmost = (np.abs(CpVals-320e-12)).argmin()
#CpVals = CpVals[leftmost:rightmost]

#leftmost = (np.abs(CsVals-55e-12)).argmin()
#rightmost = (np.abs(CsVals-75e-12)).argmin()
#CsVals = CsVals[leftmost:rightmost]
FREQ = 25*MHz_

Gammas = []

ct = TomasCircuit()


with alive_bar(len(CaVals)*len(CsVals)*len(CpVals),title="making a Gamma(Ca,Cs,Cp) movie ",length=20,bar="filling",spinner="dots_waves2") as bar:
    for i,CaVal in enumerate(CaVals):
        for CsVal in CsVals:
            for CpVal in CpVals:
                ct.set('Ca.Cs',CaVal)
                ct.set('Cs.Cs',CsVal)
                ct.set('Cp.Cs',CpVal)

                Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V3toV2'])

                Vf = abs(Solution['V3toV2.V3'][2])
                Vr = abs(Solution['V3toV2.V3'][3])

                Gammas.append(Vr/Vf)
                bar()

        X, Y = np.meshgrid(CpVals/pF_,CsVals/pF_)
        Z = np.array(Gammas).reshape(len(CsVals),len(CpVals))
        levels = np.linspace(0.0, 1.0, 11)
        CS = plt.contourf(X, Y, Z, levels=levels, cmap=cm.coolwarm, extend='min')
        colorbar = plt.colorbar(CS)
        plt.xlabel("Cp (pF)")
        plt.ylabel("Cs (pF)")
        plt.title(f"Ca = {round(CaVal/pF_,2)}")
        #plt.savefig(f"VideoRegionScan/Iframe-{i}.png")
        #plt.savefig(f"VideoRegionScan/ICa_{int(CaVal/pF_)}pF.pdf")
        plt.show()
        plt.clf()
        Gammas = []
