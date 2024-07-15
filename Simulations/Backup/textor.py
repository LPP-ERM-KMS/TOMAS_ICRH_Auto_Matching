import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs

def CalcCa(f):
    return -3.3094e-6*f**5 + 5.17212e-4*f**4-5.10998e-2*f**3+3.5211*f**2-132.074*f+2000.87

def findClosestValue(givenList, target):
    def difference(givenList):
        return abs(givenList - target)
    result = min(givenList, key=difference)
    return result

CsVals = np.linspace(25e-12,1000e-12,100000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,100000)#capacitor is defined in steps

FREQ = 25*10**6

CaVal = CalcCa(FREQ*(10**(-6)))*(10**(-12))

#initial values:
CsVal = CsVals[10]
CpVal = CpVals[10]
CsFactor = 0.2e-17
CpFactor = 0.2e-17

PlotCs = []
PlotCp = []
admittances = []

PlotCs.append(CsVal)
PlotCp.append(CpVal)

PlotEpsB = []
PlotEpsG = []

Steps = []
step = 1
Steps.append(step)

#---------#
# Circuit #
#---------#
TRLPart1 = rfTRL(L=0.235, OD=0.041, ID=0.017, dx=500)
TRLPart2 = rfTRL(L=2, OD=0.041, ID=0.017, dx=500)
TRLPart3 = rfTRL(L=3.5, OD=0.041, ID=0.017, dx=500)
#load S matrix of antenna
Antenna = rfObject(touchstone='tomas_icrh_linear_2017-vacuum.s2p')
RLCLeft = rfRLC(Cs=CsVal,Cp=CpVal)
CaRight = rfRLC(Cs=CaVal)

ct = rfCircuit(Id="MatchingCircuit")

ct.addblock('TL1',TRLPart1,ports=['A','B'],relpos=0)
ct.addblock('TL2',TRLPart2,ports=['B','C'],relpos=TRLPart1.L)
ct.addblock('TL3',TRLPart3,ports=['C','D'],relpos=TRLPart1.L + TRLPart2.L)
ct.addblock('RLoss',rfRLC(Rs=1),ports=['G','H'],relpos=TRLPart1.L + TRLPart2.L)

EndOfLine =TRLPart1.L + TRLPart2.L + TRLPart3.L 

ct.addblock('Matching',RLCLeft,ports=['F','G'],relpos=EndOfLine)
ct.addblock('Antenna',Antenna,ports=['G','H'],relpos=1+EndOfLine)
ct.addblock('PreMatch',CaRight,ports=['H','I'],relpos=2+EndOfLine)

ct.connect('TL1.A','Source')
ct.connect('TL1.B','TL2.B')
ct.connect('TL2.C','TL3.C')
ct.connect('TL3.D','Matching.F')
ct.connect('Matching.G','RLoss.G')
ct.connect('RLoss.H','Antenna.G')
ct.connect('Antenna.H','PreMatch.H')
ct.terminate('PreMatch.I',Z=0) #Grounding

Converged = False

while not Converged:
    ct.set('Matching.Cs',CsVal)
    ct.set('Matching.Cp',CpVal)

    Solution = ct.Solution(f=FREQ, E={'Source': 360000}, nodes=['TL1.A','TL2.B','TL3.C','PreMatch.D'])

    """
    Calculate new values of Cs and Cp
    """
    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([0.5,2])
    voltages = np.zeros(len(MeasurePoints))
    Vs = np.zeros(len(MeasurePoints))

    #constants:
    beta = 2*np.pi*FREQ/(3*(10**8))
    S = np.sin(2*beta*MeasurePoints) #array
    C = np.cos(2*beta*MeasurePoints) #array

    #Make array of voltages on the various points:
    voltages[0] = np.abs(Solution['TL2.B'][0])
    voltages[1] = np.abs(Solution['TL3.C'][0])

    Vf = abs(Solution['TL1.A'][2])
    Vr = abs(Solution['TL1.A'][3])
    print("reflection coeff:")
    print(np.abs(Vr/Vf))
    ReflCoeff = Vr/Vf

    EpsA = (voltages[0] - Vf)*S[1] - (voltages[1] - Vf)*S[0]
    EpsG = (voltages[0] - Vf)*C[1] - (voltages[1] - Vf)*C[0]
    print("Eps:")
    print(EpsA + 1j*EpsG)

    admittance = EpsG+1 + 1j*EpsA
    admittances.append(admittance)
    
    DeltaCs = CsFactor*EpsG
    DeltaCp = CpFactor*EpsA

    CsVal += DeltaCs
    CpVal += DeltaCp
        
    IndexCloseToNewCs = np.where(findClosestValue(CsVals, CsVal) == CsVals)#get index of closest lying computed 
    IndexCloseToNewCp = np.where(findClosestValue(CpVals, CpVal) == CpVals)#get index of closest lying computed 

    CsVal = CsVals[IndexCloseToNewCs][0]
    CpVal = CpVals[IndexCloseToNewCp][0]

    PlotCs.append(CsVal)
    PlotCp.append(CpVal)

    PlotEpsB.append(EpsA)
    PlotEpsG.append(EpsG)

    #Converged enough:
    if (np.abs(ReflCoeff) < 0.01):
        print("nicely converged")
        maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': 360000})
        Converged = True
    step += 1
    Steps.append(step)
    if ((PlotCs[-2] == PlotCs[-1]) and (PlotCp[-2] == PlotCp[-1])):
        print("converged due to a lack of changes")
        maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': 360000})
        Converged = True

plotVSWs(VSWs) 
plt.show()

print("final Gamma:")
print(np.abs(ReflCoeff))
plt.plot(np.arange(0,len(admittances)),np.real(admittances),label="real part")
plt.plot(np.arange(0,len(admittances)),np.imag(admittances),label="imaginary part")
plt.show()
plt.plot(Steps,PlotCs,label="Cs")
plt.plot(Steps,PlotCp,label="Cp")
plt.xlabel("step")
plt.ylabel("Capacitor value (F)")
plt.legend()
plt.show()

plt.plot(Steps[1:],PlotEpsB,label="$\epsilon_B$")
plt.plot(Steps[1:],PlotEpsG,label="$\epsilon_G$")
plt.ylabel("$\epsilon$")
plt.xlabel("step")
plt.legend()
plt.show()
