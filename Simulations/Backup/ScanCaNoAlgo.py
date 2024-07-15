import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import Circuit, Circuit2
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
#This will take too long, so we'll limit the region
leftmost = (np.abs(CaVals-100e-12)).argmin()
rightmost = (np.abs(CaVals-300e-12)).argmin()
CaVals = CaVals[leftmost:rightmost]


FREQ = 25.0*MHz_
Freqs = np.linspace(25*MHz_,49*MHz_,25)

Factor = 100
Epss = [0,1,2,3,4]
EpsLabels = ["eps=100","eps=200","eps=500","eps=842","eps=1000","eps=1500","eps=2000","eps=2500","eps=3000"]

for FREQ in Freqs:
    print(f"Matching at {FREQ/MHz_}MHz")
    for eps in Epss:
        ct = Circuit2(WhichStrap=eps)
        with alive_bar(len(CaVals),title=f"Matching {EpsLabels[eps]}",length=20,bar="filling",spinner="dots_waves2") as bar:
            Gammas = []
            CsFactor = Factor*pF_
            CpFactor = Factor*pF_
            for CaVal in CaVals:
                Steps = [1]
                #ct.set('PlasmaLike.Rs',Plasma)
                ct.set('Ca.Cs',CaVal)
                Converged = False
                #initial values:
                CsVal = 750*pF_ #CsVals[100]
                CpVal = 750*pF_ #CpVals[100]

                IntermediaryResult = []


                while not Converged:
                    ct.set('Cs.Cs',CsVal)
                    ct.set('Cp.Cs',CpVal)

                    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['CptoV0','V0toV1','V2toV3'])

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

                    Gamma = np.abs(Vr/Vf)


                    Vs = (np.abs(V)**2)/(np.abs(Vf)**2) 

                    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
                    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

                    EpsB = 2*v/((1+u)**2 + v**2)
                    EpsG = (1 - ((1-u**2-v**2)/((1+u)**2 + v**2)))

                    yCp = -Solution['CptoV0.Cp'][1]/Solution['CptoV0.Cp'][0]*50
                    dCmax = 10e-12
                    dCs = max(-dCmax, min((yCp.real-1)*1, dCmax))
                    dCp = max(-dCmax, min(yCp.imag*1, dCmax))

                    CsVal += dCs
                    CpVal += dCp
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
                    
                    """
                    if Gamma < 0.1:
                        print(f"freq {FREQ/1e6}MHz is converging!")
                        IntermediaryResult.append([Gamma,Factor,CaVal,CsVal,CpVal])

                    """

                bar()
                IntermediaryResult = np.array(IntermediaryResult)
                Gammas.append(IntermediaryResult[np.where(IntermediaryResult[:,0] == min(IntermediaryResult[:,0])),0][0][0])

        plt.xlabel("Ca (pF)")
        plt.ylabel("$1 - |\Gamma$|")
        plt.plot(CaVals/pF_,1-np.array(Gammas),label=EpsLabels[eps])
        plt.title(f"Scan of Ca Values at {FREQ/MHz_}MHz, each time matching using the 4V matching algorithm")
    plt.legend()
    plt.show()
