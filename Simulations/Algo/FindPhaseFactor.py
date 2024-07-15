import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod as TomasCircuit
from alive_progress import alive_bar
from Algos import *
from smithplot import SmithAxes


MHz_ = 1e6
pF_ = 1e-12

phasefactors = np.linspace(-2,-0.5,4)
resolution = 10
maxcs = 1000e-12
mincs = 7e-12
mincp = 35e-12
maxcp = 1000e-12

FREQ = 25*MHz_


CsVals = np.linspace(mincs,maxcs,resolution)#capacitor is defined in steps
StartCsVals = np.linspace(200e-12,800e-12,10)#capacitor is defined in steps
#leftmost = (np.abs(CsVals-120e-12)).argmin()
#rightmost = (np.abs(CsVals-170e-12)).argmin()
#CsVals = CsVals[leftmost:rightmost]

CpVals = np.linspace(mincp,maxcp,resolution)#capacitor is defined in steps
StartCpVals = np.linspace(200e-12,800e-12,10)#capacitor is defined in steps
#leftmost = (np.abs(CpVals-100e-12)).argmin()
#rightmost = (np.abs(CpVals-750e-12)).argmin()
#CpVals = CpVals[leftmost:rightmost]

U = []
V = []

def MatchAble(algorithm,ct,CpVal,CsVal,FREQ,ProbeIndexA=None,ProbeIndexB=None,ProbeIndexC=None,phasefactor=None):
    i = 0
    matched = False
    while ((i < 1500) and not matched):
        
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

        if abs(Vr/Vf) <= 0.1:
            matched = True
        if phasefactor:
            EpsG,EpsB = algorithm(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC,phasefactor)
        else:
            EpsG,EpsB = algorithm(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC)
        
        CsVal += EpsG*1*pF_
        CpVal += EpsB*1*pF_

        i += 1

    return matched


Gammas = []
ratios = []

ct = TomasCircuit()

CaVal = 200*pF_
ct.set('Ca.Cs',CaVal)
num = 0
tot = 0
avg_steps = 0

with alive_bar(len(CsVals)*len(CpVals)*len(phasefactors),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for phasefactor in phasefactors:
        for j,CsVal in enumerate(StartCsVals):
            for k,CpVal in enumerate(StartCpVals):
                ct.set('Cs.Cs',CsVal)
                ct.set('Cp.Cs',CpVal)
                Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V3toV2'])
                steps,matched = MatchAble(ModAlgo3V,ct,CpVal,CsVal,FREQ,0,2,3,phasefactor)
                tot+=1
                avg_steps += steps/(len(StartCsVals)*len(StartCpVals))
                if matched:
                    num+=1
                bar()
        ratios.append(num/tot)
        num =0
        tot =0

plt.plot(phasefactors,ratios)
plt.xlabel("phasefactor")
plt.ylabel("ratios")
plt.show()
