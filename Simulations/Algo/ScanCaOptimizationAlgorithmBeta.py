import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurCircuit as TomasCircuit
from alive_progress import alive_bar
from Algos import *

MHz_ = 1e6
pF_ = 1e-12

global mincs = 25e-12
global maxcs = 1000e-12
global mincp = 7e-12
global maxcp = 1000e-12

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
def Match(algorithm,ct,CpVal,CsVal,FREQ):
    matched = False
    CsVals = np.linspace(mincs,maxcs,5000)#capacitor is defined in steps
    CpVals = np.linspace(mincp,maxcp,5000)#capacitor is defined in steps
    j=0
    while not matched and j < 100:
        CsVal = np.random.choice(CsVals)
        CpVal = np.random.choice(CpVals)
        while ((i < 10000) and not matched):
            
            if CsVal > maxcs or CsVal < mincs:
                break
            if CpVal > maxcp or CpVal < mincp:
                break

            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)

            Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V1toV0','V3toV2'])

            Vf = Solution['V3toV2.V3'][2]
            Vr = Solution['V3toV2.V3'][3]
            Gamma = Vr/Vf

            if abs(Vr)/abs(Vf) <= 0.05:
                matched = True
                return CsVal,CpVal

            EpsG,EpsB = algorithm(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC)
            
            CsVal += EpsG*1*pF_
            CpVal += EpsB*1*pF_

            i += 1
        j+=1
    return 100e-12,100e-12



FREQ = 25.0*MHz_
Gammas = []
ct = TomasCircuit()
CaVals = np.linspace(250e-12,350e-12,50)#capacitor is defined in steps
FCsVals = []
FCpVals = []
CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps

with alive_bar(len(CaVals),title="Scanning Ca and matching",length=20,bar="filling",spinner="dots_waves2") as bar:
    for CaVal in CaVals:
        CsVal,CpVal = Match(DirCoupler,ct,CpVal,CsVal,FREQ):
        result = optimize.minimize(ComputeGamma, np.array([CsVal,CpVal]),args=(CaVal,FREQ),bounds=bounds,method='Nelder-Mead')
        Gamma = result.fun
        CsVal = result.x[1]
        CpVal = result.x[0]

        Gammas.append(Gammas)
        FCsVals.append(CsVal)
        FCpVals.append(CpVal)
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
