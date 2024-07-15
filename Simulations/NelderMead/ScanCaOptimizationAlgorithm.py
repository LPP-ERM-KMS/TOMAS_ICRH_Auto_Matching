import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurCircuit as TomasCircuit
from alive_progress import alive_bar

MHz_ = 1e6
pF_ = 1e-12

def findClosestValue(givenList, target):
    def difference(givenList):
        return abs(givenList - target)
    result = min(givenList, key=difference)
    return result

def ComputeGamma(CVals,CaVal,FREQ=25*MHz_):
    CpVal = CVals[0]
    CsVal = CVals[1]

    ct.set('PreMatch.Cs',CaVal)
    ct.set('Matching.Cs',CsVal)
    ct.set('Matching.Cp',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V3toV2'])

    Vf = abs(Solution['V3toV2.V3'][2])
    Vr = abs(Solution['V3toV2.V3'][3])

    Gamma = Vr/Vf
    return Gamma



FREQ = 35.0*MHz_
Gammas = []
ct = TomasCircuit()
CaVals = np.linspace(250e-12,350e-12,50)#capacitor is defined in steps
FCsVals = []
FCpVals = []
CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps

with alive_bar(len(CaVals),title="Scanning Ca and matching",length=20,bar="filling",spinner="dots_waves2") as bar:
    for CaVal in CaVals:
        localGammas = np.zeros(10)
        localFCsVals= np.zeros(10)
        localFCpVals= np.zeros(10)
        bounds = optimize.Bounds(lb=np.array([25e-12,7e-12]),ub=np.array([1000e-12,1000e-12]))
        for i in range(10):
            result = optimize.minimize(ComputeGamma, np.array([np.random.choice(CsVals),np.random.choice(CpVals)]),args=(CaVal,FREQ),bounds=bounds,method='Nelder-Mead')
            localGammas[i] = result.fun
            localFCsVals[i] = result.x[1]
            localFCpVals[i] = result.x[0]

        Gammas.append(min(localGammas))
        index = np.where(min(localGammas) == localGammas)
        FCsVals.append(localFCsVals[index[0][0]])
        FCpVals.append(localFCpVals[index[0][0]])
        bar()

plt.plot(CaVals/pF_,Gammas)
plt.title("Ca scan, matching using Nelder-Mead")
plt.xlabel("Ca (pF)")
plt.ylabel("$\Gamma$")
plt.show()

print(FCpVals)
plt.plot(CaVals/pF_,np.array(FCpVals)/pF_,label="Cp")
plt.title("Ca scan, matching using Nelder-Mead")
plt.xlabel("Ca (pF)")
plt.ylabel("C (pF)")
plt.plot(CaVals/pF_,np.array(FCsVals)/pF_,label="Cs")
plt.legend()
plt.show()
