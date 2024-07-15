import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import Circuit
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

ct = Circuit()

CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps
CaVals = np.linspace(35e-12,1000e-12,8000)#capacitor is defined in steps
#This will take too long, so we'll limit the region
leftmost = (np.abs(CaVals-150e-12)).argmin()
rightmost = (np.abs(CaVals-180e-12)).argmin()
CaVals = CaVals[leftmost:rightmost]

#file = open('Factors.csv','w')
file = open('Test.csv','w')
writer = csv.writer(file)

writer.writerow(["FREQ","Factor","CaVal","CsVal","CpVal"])

FREQ = 25*MHz_
#frequencies = np.linspace(25,55,30)*MHz_
frequencies = [50*MHz_]

Factors = np.linspace(25,225,50)
StepsFactors = np.zeros(len(Factors))

with alive_bar(len(CaVals)*len(Factors)*len(frequencies),title="Matching",length=20,bar="filling",spinner="dots_waves2") as bar:
    for FREQ in frequencies:
        for FactorIter,Factor in enumerate(Factors):
            FinalGammas = []
            CsFactor = Factor*pF_
            CpFactor = Factor*pF_
            for CaVal in CaVals:
                Steps = [1]
                #ct.set('PlasmaLike.Rs',Plasma)
                ct.set('Ca.Cs',CaVal)
                Converged = False
                #initial values:
                CsVal = 133*pF_ #CsVals[100]
                CpVal = 47.94*pF_ #CpVals[100]

                Gammas = []
                while not Converged:
                    ct.set('Cs.Cs',CsVal)
                    ct.set('Cp.Cs',CpVal)

                    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

                    ######################################
                    # Calculate new values of Cs and Cp  #
                    ######################################

                    # Get voltages at positions in Measurepoints:
                    MeasurePoints = np.array([0.235,0.895,1.69,2.35])
                    V = np.zeros(len(MeasurePoints))

                    #constants:
                    beta = 2*np.pi*FREQ/(3*(10**8))
                    S = np.sin(2*beta*MeasurePoints) #array
                    C = np.cos(2*beta*MeasurePoints) #array
                    BigS = (S[0] - S[1])/(S[2] - S[3])
                    BigC = (C[0] - C[1])/(C[2] - C[3])

                    #Make array of voltages on the various points:
                    V[0] = np.abs(Solution['V0toV1.V0'][0])
                    V[1] = np.abs(Solution['V0toV1.V1'][0])
                    V[2] = np.abs(Solution['V2toV3.V2'][0])
                    V[3] = np.abs(Solution['V2toV3.V3'][0])

                    Vf = abs(Solution['V2toV3.V3'][2])
                    Vr = abs(Solution['V2toV3.V3'][3])

                    Gammas.append(Vr/Vf)

                    Vs = (np.abs(V)**2)/(np.abs(Vf)**2) 

                    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
                    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

                    EpsB = 2*v/((1+u)**2 + v**2)
                    EpsG = (1 - ((1-u**2-v**2)/((1+u)**2 + v**2)))

                    CsVal += CsFactor*EpsG
                    CpVal += CpFactor*EpsB
                    CpVal = CpVals[(np.abs(CpVals-CpVal)).argmin()]
                    CsVal = CsVals[(np.abs(CsVals-CsVal)).argmin()]

                    #Converged enough:
                    if (np.abs(Gammas[-1]) < 0.02):
                        print("nicely converged")
                        StepsFactors[FactorIter] = Steps[-1]
                        writer.writerow([FREQ,Factor,CaVal,CsVal,CpVal])
                        Converged = True

                    Steps.append(Steps[-1]+1)
                    if Steps[-1] > 10:
                        Converged = True
                        if StepsFactors[FactorIter] == 0:
                            StepsFactors[FactorIter] = None
                bar()
file.close()

"""
plt.legend()
plt.xlabel("Ca (pF)")
plt.ylabel("$|\Gamma$|")
plt.show()

print(StepsFactors)
plt.scatter(Factors,StepsFactors)
plt.title("25 MHz")
plt.xlabel("$C=S$")
plt.ylabel("Steps")
plt.show()
"""
