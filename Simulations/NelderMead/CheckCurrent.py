import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from scipy.stats import qmc
from TomasCircuit import ArthurMod as TomasCircuit
from alive_progress import alive_bar

MHz_ = 1e6
pF_ = 1e-12

CsVals = np.linspace(25e-12,1000e-12,1000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,1000)#capacitor is defined in steps

leftmost = (np.abs(CpVals-800e-12)).argmin()
rightmost = (np.abs(CpVals-1000e-12)).argmin()
CpVals = CpVals[leftmost:rightmost]

leftmost = (np.abs(CsVals-800e-12)).argmin()
rightmost = (np.abs(CsVals-1000e-12)).argmin()
CsVals = CsVals[leftmost:rightmost]


FREQ = 25*MHz_

Gammas = []
Iratios = []

ct = TomasCircuit()

CaVal = 275.72*pF_

with alive_bar(len(CsVals)*len(CpVals),title="making a Gamma(Cs,Cp) plot and a max(ICa(Cp,Cs),ICp(Cp,Cs),ICs(Cp,Cs)) plot",length=20,bar="filling",spinner="dots_waves2") as bar:
#for CaVal in CaVals:
    for CsVal in CsVals:
        for CpVal in CpVals:
            ct.set('Ca.Cs',CaVal)
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)

            Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*1.5*1E3)}, nodes=['V1toV0','V3toV2','Cs','Cp','Ca'])
            # power range: 1000-6000W i.e 10log_10(1000/10e-3W)dBm - 10log_10(6000/10e-3W)dBm
            # Or 115.13dBm - 133.0468dBm
            # But pyRFtk works with E=Vpeak, i.e E=Vpeak= sqrt(2*P_avg*R) = sqrt(2*Z_char*P_avg)

            rICa = (1/np.sqrt(2))*np.abs(Solution['Ca.cap'][1])/138
            rICp = (1/np.sqrt(2))*np.abs(Solution['Cp.Cp'][1])/70
            rICs = (1/np.sqrt(2))*np.abs(Solution['Cs.cap'][1])/140

            Iratios.append(max(np.array([rICa,rICp,rICs])))

            Vf = abs(Solution['V3toV2.V3'][2])
            Vr = abs(Solution['V3toV2.V3'][3])

            Gammas.append(Vr/Vf)

            if Vr/Vf < 0.1:
                print(f"Cs:{CsVal} Cp:{CpVal} rICa:{rICa} rICp:{rICp} rICs:{rICs}")
            bar()

    X, Y = np.meshgrid(CpVals/pF_,CsVals/pF_)
    Z = np.array(Gammas).reshape(len(CsVals),len(CpVals))
    levels = np.linspace(0.0, 1.0, 11)
    CS = plt.contourf(X, Y, Z, levels=levels, cmap=cm.coolwarm, extend='min')
    colorbar = plt.colorbar(CS)
    plt.xlabel("Cp (pF)")
    plt.ylabel("Cs (pF)")
    plt.title(f"Ca = {round(CaVal/pF_,2)}")
    plt.show()

    X, Y = np.meshgrid(CpVals/pF_,CsVals/pF_)
    Z = np.array(Iratios).reshape(len(CsVals),len(CpVals))
    levels = np.linspace(0.0, 1.1, 12)
    CS = plt.contourf(X, Y, Z, levels=levels, cmap=cm.coolwarm, extend='max')
    colorbar = plt.colorbar(CS)
    plt.xlabel("Cp (pF)")
    plt.ylabel("Cs (pF)")
    plt.title(f"Ca = {round(CaVal/pF_,2)}")
    plt.show()
    Iratios = []
