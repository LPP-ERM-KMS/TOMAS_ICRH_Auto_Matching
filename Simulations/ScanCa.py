import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod
from alive_progress import alive_bar
from Algos import *

MHz_ = 1e6
pF_ = 1e-12

def findClosestValue(givenList, target):
    def difference(givenList):
        return abs(givenList - target)
    result = min(givenList, key=difference)
    return result

def GetCa(f):
    return np.exp(-0.1004086*f-20.0279) + 2.26764*(10**(-11))


CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps
CaVals = np.linspace(35e-12,1000e-12,800)#capacitor is defined in steps
#This will take too long, so we'll limit the region
#leftmost = (np.abs(CaVals-35e-12)).argmin()
#rightmost = (np.abs(CaVals-1000e-12)).argmin()
#CaVals = CaVals[leftmost:rightmost]


Freqs = np.linspace(25*MHz_,49*MHz_,2)

Factor = 10

with alive_bar(len(CaVals)*len(Freqs),title="Matching",length=20,bar="filling",spinner="dots_waves2") as bar:
    for FREQ in Freqs:
        print(f"Matching {FREQ/MHz_}MHz")

        #ct = TomasCircuit()
        ct = ArthurMod()
        Gammas = []
        CsFactor = Factor*pF_
        CpFactor = Factor*pF_
        for CaVal in CaVals:
            Steps = [1]
            #ct.set('PlasmaLike.Rs',Plasma)
            ct.set('Ca.Cs',CaVal)
            Converged = False
            #initial values:
            CsVal = 100*pF_ #CsVals[100]
            CpVal = 100*pF_ #CpVals[100]

            IntermediaryResult = []


            while not Converged:
                ct.set('Cs.Cs',CsVal)
                ct.set('Cp.Cs',CpVal)

                Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

                Vf = np.abs(Solution['V3toV2.V3'][2])
                Vr = np.abs(Solution['V3toV2.V3'][3])

                Gamma = Vr/Vf

                EpsG,EpsB = ModAlgo3V(Solution,FREQ,0,1,3)

                CsVal += CsFactor*EpsG
                CpVal += CpFactor*EpsB
                CpVal = CpVals[(np.abs(CpVals-CpVal)).argmin()]
                CsVal = CsVals[(np.abs(CsVals-CsVal)).argmin()]

                #Converged enough:
                if Gamma < 0.1:
                    print("nicely converged")
                    Converged = True
                    IntermediaryResult.append([Gamma,Factor,CaVal,CsVal,CpVal])
                
                elif Steps[-1] > 100:
                    Converged = True
                    IntermediaryResult.append([Gamma,Factor,CaVal,CsVal,CpVal])

                Steps.append(Steps[-1]+1)
                
            bar()
            IntermediaryResult = np.array(IntermediaryResult)
            Gammas.append(IntermediaryResult[np.where(IntermediaryResult[:,0] == min(IntermediaryResult[:,0])),0][0][0])

        plt.xlabel("Ca (pF)")
        plt.ylabel("$|\Gamma$|")
        plt.plot(CaVals/pF_,np.array(Gammas))

        plt.title(f"Scan of Ca Values at {FREQ/MHz_}MHz, each time matching using the 3V matching algorithm")
        plt.legend()
        plt.show()
