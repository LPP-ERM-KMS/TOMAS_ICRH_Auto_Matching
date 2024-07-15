import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import SimpleCircuit
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


CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps
CaVals = np.linspace(35e-12,1000e-12,8000)#capacitor is defined in steps

FREQ = 25.0*MHz_
#Freqs = np.linspace(11*MHz_,49*MHz_,39)
Freqs = np.linspace(25*MHz_,49*MHz_,2)

Factor = 100

with alive_bar(len(CaVals)*len(Freqs),title="Matching",length=20,bar="filling",spinner="dots_waves2") as bar:
    for FREQ in Freqs:
        print(f"Matching {FREQ/MHz_}MHz")
        ct = SimpleCircuit()
        Gammas = []
        CsFactor = Factor*pF_
        CpFactor = Factor*pF_
        for CaVal in CaVals:
            Steps = [1]
            #ct.set('PlasmaLike.Rs',Plasma)
            ct.set('PreMatch.Cs',CaVal)
            Converged = False
            #initial values:
            CsVal = 100*pF_ #CsVals[100]
            CpVal = 100*pF_ #CpVals[100]

            IntermediaryResult = []


            while not Converged:
                ct.set('Matching.Cs',CsVal)
                ct.set('Matching.Cp',CpVal)

                Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V1toV0','V3toV2'])

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
                V[0] = np.abs(Solution['V1toV0.V0'][0])
                V[1] = np.abs(Solution['V1toV0.V1'][0])
                V[2] = np.abs(Solution['V3toV2.V2'][0])
                V[3] = np.abs(Solution['V3toV2.V3'][0])

                Vf = abs(Solution['V3toV2.V3'][2])
                Vr = abs(Solution['V3toV2.V3'][3])

                Gamma = np.abs(Vr/Vf)


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
                if Gamma < 0.1:
                    print("nicely converged")
                    Converged = True
                    IntermediaryResult.append([Gamma,Factor,CaVal,CsVal,CpVal])
                
                elif Steps[-1] > 20:
                    Converged = True
                    IntermediaryResult.append([Gamma,Factor,CaVal,CsVal,CpVal])

                Steps.append(Steps[-1]+1)
                
            bar()
            IntermediaryResult = np.array(IntermediaryResult)
            Gammas.append(IntermediaryResult[np.where(IntermediaryResult[:,0] == min(IntermediaryResult[:,0])),0][0][0])

        plt.xlabel("Ca (pF)")
        plt.ylabel("$1 - |\Gamma$|")
        plt.plot(CaVals/pF_,1-np.array(Gammas))
        plt.title(f"Scan of Ca Values at {FREQ/MHz_}MHz, each time matching using the 4V matching algorithm")
        plt.legend()
        plt.show()
        #plt.savefig(f"figures/Ca_Scan_{FREQ/MHz_}_MHz.pdf")
        #plt.figure()
