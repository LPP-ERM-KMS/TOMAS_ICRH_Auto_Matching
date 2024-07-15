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

ct = Circuit()

CsVals = np.linspace(25e-12,1000e-12,1000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,1000)#capacitor is defined in steps

MHz_ = 1e6
pF_ = 1e-12

FREQ = 40*MHz_

CaVals=np.linspace(120,220,1000)*pF_ #CalcCa(FREQ*(10**(-6)))*(10**(-12))

CsFactor = 0.01*pF_
CpFactor = 0.01*pF_

PlotEpsB = []
PlotEpsG = []

Plasmas = np.linspace(0,1,5)

for Plasma in Plasmas:
    Steps = []
    step = 1
    Steps.append(step)
    ct.set('PlasmaLike.Rs',Plasma)
    SpecialRs = []
    SpecialCas = []
    admittances = []
    ReflCoeffs = []
    for CaVal in CaVals:
        ct.set('Ca.Cs',CaVal)
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

            #ct = Circuit(CaVal,CpVal,CsVal)
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

            EpsB = (V[1] - Vf)*S[3] - (V[3] - Vf)*S[1]
            EpsG = (V[1] - Vf)*C[3] - (V[3] - Vf)*C[1]

            print("Eps:")
            print(EpsB + 1j*EpsG)
            admittance = EpsG+1 + 1j*EpsB
            admittances.append(admittance)
            
            #DeltaCs = CsFactor*EpsG
            #DeltaCp = CpFactor*EpsB
            gm = ct.getS(25e6)[0]
            tSol = ct.Solution(25e6, {'Source':np.sqrt(2*50*5*1E3)}, nodes=['CptoV0.Cp'])
            yCp = -tSol['CptoV0.Cp'][1]/tSol['CptoV0.Cp'][0]*50
            DeltaCs = max(-10E-12, min((yCp.real-1), 10E-12))
            DeltaCp = max(-10E-12, min(yCp.imag*(-1), 10E-12))


            CsVal += DeltaCs
            CpVal += DeltaCp

            if CsVal > 1000e-12:
                CsVal = 1000e-12
            if CsVal < 7e-12:
                CsVal = 7e-12
            if CpVal > 1000e-12:
                CpVal = 1000e-12
            if CpVal < 7e-12:
                CpVal = 7e-12

                
            #IndexCloseToNewCs = np.where(findClosestValue(CsVals, CsVal) == CsVals)#get index of closest lying computed 
            #IndexCloseToNewCp = np.where(findClosestValue(CpVals, CpVal) == CpVals)#get index of closest lying computed 

            #CsVal = CsVals[IndexCloseToNewCs][0]
            #CpVal = CpVals[IndexCloseToNewCp][0]

            PlotCs.append(CsVal)
            PlotCp.append(CpVal)

            PlotEpsB.append(EpsB)
            PlotEpsG.append(EpsG)

            #Converged enough:
            if (np.abs(ReflCoeff) < 0.01):
                print("nicely converged")
                Converged = True
                SpecialRs.append(ReflCoeff)
                SpecialCas.append(CaVal)
            step += 1
            Steps.append(step)
            if ((PlotCs[-2] == PlotCs[-1]) and (PlotCp[-2] == PlotCp[-1])):
                print("converged due to a lack of changes")
                Converged = True
                CpFactor *= -0.8
            if step > 500:
                maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': np.sqrt(2*50*6*1E3)})
                Converged = True
        ReflCoeffs.append(ReflCoeff)


    #plt.ylim(0,0.5)
    plt.plot(np.array(CaVals)/pF_,ReflCoeffs,label=f"Proxy For Plasma Density: Rs={Plasma}")
    plt.plot(np.array(SpecialCas)/pF_,SpecialRs,color="red")
plt.legend()
plt.xlabel(r"$C_a (pF)$")
plt.ylabel(r"$\Gamma$")
plt.title(f"Pre-matching capacitor scan at {FREQ/MHz_}MHz")
plt.xlim(125,220)
plt.show()

