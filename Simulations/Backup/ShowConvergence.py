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


ct = Circuit2(WhichStrap=5)

CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps
CaVals = np.linspace(35e-12,1000e-12,8000)#capacitor is defined in steps
#This will take too long, so we'll limit the region
leftmost = (np.abs(CaVals-35e-12)).argmin()
rightmost = (np.abs(CaVals-170e-12)).argmin()
CaVals = CaVals[leftmost:rightmost]

FREQ = 11.0*MHz_

Factor = 30

Gammas = []
CsFactor = Factor*pF_
CpFactor = Factor*pF_
IntermediaryResult = []

CaVal = 982.9*pF_
ct.set('Ca.Cs',CaVal)
Converged = False
#initial values:
CsVal = 1000*pF_ #CsVals[100]
CpVal = 1000*pF_ #CpVals[100]


PlotEpsB = []
PlotEpsG = []
PlotCs = []
PlotCp = []
PlotCpI = []
PlotCsI = []
PlotCaI = []
Steps = []
i = 0

while not Converged:
    Steps.append(i)
    i += 1
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)

    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3','Cp','Cs','Ca'])

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

    CpI = np.abs(Solution['Cp.sc'][1])
    CsI = np.abs(Solution['Cs.cap'][1])
    CaI = np.abs(Solution['Ca.cap'][1])

    Vf = abs(Solution['V2toV3.V3'][2])
    Vr = abs(Solution['V2toV3.V3'][3])

    Gammas.append(Vr/Vf)

    #Converged enough:
    if Gammas[-1] < 0.1:
        print("nicely converged")
        Converged = True
    
    elif i > 20:
        Converged = True


    PlotCs.append(CsVal/pF_)
    PlotCp.append(CpVal/pF_)


    Vs = (np.abs(V)**2)/(np.abs(Vf)**2) 

    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

    EpsB = 2*v/((1+u)**2 + v**2)#*Gammas[-1]
    EpsG = (1 - ((1-u**2-v**2)/((1+u)**2 + v**2)))#*Gammas[-1]

    
    CsVal += CsFactor*EpsG
    CpVal += CpFactor*EpsB
    CpVal = CpVals[(np.abs(CpVals-CpVal)).argmin()]
    CsVal = CsVals[(np.abs(CsVals-CsVal)).argmin()]

    PlotEpsB.append(EpsB)
    PlotEpsG.append(EpsG)
    PlotCpI.append(CpI)
    PlotCsI.append(CsI)
    PlotCaI.append(CaI)

plt.xlabel("Steps")
plt.ylabel("$|\Gamma$|")
plt.plot(Steps,Gammas)
plt.show()

plt.xlabel("Steps")
plt.ylabel("$|\epsilon$|")
plt.plot(Steps,PlotEpsB,label="EpsB ($\propto$ Cp)")
plt.plot(Steps,PlotEpsG,label="EpsG ($\propto$ Cs)")
plt.legend()
plt.show()

plt.xlabel("$|\Gamma$|")
plt.ylabel("$|\epsilon$|")
plt.plot(Gammas,PlotEpsB,label="EpsB ($\propto$ Cp)")
plt.plot(Gammas,PlotEpsG,label="EpsG ($\propto$ Cs)")
plt.legend()
plt.show()

plt.xlabel("Steps")
plt.ylabel("C (pF)")
plt.plot(Steps,PlotCs,label="Cs")
plt.plot(Steps,PlotCp,label="Cp")
plt.legend()
plt.show()

plt.xlabel("Steps")
plt.ylabel("Current as % of maximum")
plt.plot(Steps,np.array(PlotCpI)/70.0,label="Cp")
plt.plot(Steps,np.array(PlotCsI)/140.0,label="Cs")
plt.plot(Steps,np.array(PlotCaI)/138.0,label="Ca")
plt.legend()
plt.show()
