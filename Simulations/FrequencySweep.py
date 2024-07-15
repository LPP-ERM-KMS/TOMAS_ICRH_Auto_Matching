import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod as Circuit
from alive_progress import alive_bar

MHz_ = 1e6
pF_ = 1e-12

def findClosestValue(givenList, target):
    def difference(givenList):
        return abs(givenList - target)
    result = min(givenList, key=difference)
    return result

def GetCa(f):
    return np.exp(-0.1004086*f-20.0279) + 2.26764*(10**(-11))


CsVal = 200e-12
CpVal = 200e-12
CaVal = 200e-12
#This will take too long, so we'll limit the region
#leftmost = (np.abs(CaVals-35e-12)).argmin()
#rightmost = (np.abs(CaVals-1000e-12)).argmin()
#CaVals = CaVals[leftmost:rightmost]


Freqs = np.linspace(10*MHz_,55*MHz_,200)

Gammas = []
with alive_bar(len(Freqs),title="Matching",length=20,bar="filling",spinner="dots_waves2") as bar:
    for FREQ in Freqs:
        ct = Circuit()
        ct.set('Ca.Cs',CaVal)
        Converged = False
        ct.set('Cs.Cs',CsVal)
        ct.set('Cp.Cs',CpVal)
        Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V3toV2'])

        Vf = abs(Solution['V3toV2.V3'][2])
        Vr = abs(Solution['V3toV2.V3'][3])

        Gamma = np.abs(Vr/Vf)

        bar()
        Gammas.append(Gamma)

plt.xlabel("Freq (MHz)")
plt.ylabel("$\Gamma$")
plt.plot(Freqs/MHz_,np.array(Gammas))

plt.title(f"Frequency sweep for Ca={CaVal/pF_}pF, Cp={CpVal/pF_}pF and Cs={CsVal/pF_}pF")
plt.show()
