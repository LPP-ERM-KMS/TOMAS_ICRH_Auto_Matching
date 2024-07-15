import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import Circuit

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

MHz_ = 1e6
pF_ = 1e-12

FREQ = 45*MHz_

CsFactor = 10*pF_
CpFactor = 10*pF_

PlotGammaR = []
PlotGammeI = []

#Plasmas = np.linspace(0,1,5)

Steps = []
step = 1
Steps.append(step)
#ct.set('PlasmaLike.Rs',Plasma)
ct.set('Ca.Cs',GetCa(FREQ))
Converged = False
PlotCs = []
PlotCp = []
#initial values:
CsVal = 133*pF_ #CsVals[100]
CpVal = 47.94*pF_ #CpVals[100]
PlotCs.append(CsVal)
PlotCp.append(CpVal)

while not Converged:
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    ct.set('Ca.Cs',CaVal)

    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

    #Calculate new values of Cs and Cp

    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([0.235,0.66,0.66+0.795,0.235+0.66+0.795+0.66])
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

    print("reflection coeff:")
    print(np.abs(Vr/Vf))
    ReflCoeff = np.abs(Vr/Vf)

    Vs = (np.abs(V)**2)/(np.abs(Vf)**2) 

    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

    EpsB = 2*v/((1+u)**2 + v**2)
    EpsG = (1 - ((1-u**2-v**2)/((1+u)**2 + v**2)))

    CsVal += CsFactor*EpsG
    CpVal += CpFactor*EpsB
    CpVal = CpVals[(np.abs(CpVals-CpVal)).argmin()]
    CsVal = CsVals[(np.abs(CsVals-CsVal)).argmin()]
    PlotCs.append(CsVal)
    PlotCp.append(CpVal)

    #Converged enough:
    if (np.abs(ReflCoeff) < 0.01):
        print("nicely converged")
        Converged = True
    step += 1
    Steps.append(step)
    if ((PlotCs[-2] == PlotCs[-1]) and (PlotCp[-2] == PlotCp[-1])):
        print("converged due to a lack of changes")
        Converged = True
        CpFactor *= -0.8
    if step > 5000:
        maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': np.sqrt(2*50*6*1E3)})
        print("too many steps")
        Converged = True

