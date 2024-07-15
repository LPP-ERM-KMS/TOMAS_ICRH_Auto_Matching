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

CsVals = np.linspace(25e-12,1000e-12,400)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,400)#capacitor is defined in steps

FREQ = 25*10**6

CaVal = CalcCa(FREQ*(10**(-6)))*(10**(-12))

#initial values:
CsVal = CsVals[100]

PlotCs = []

PlotEpsB = []
PlotEpsG = []

#---------#
# Circuit #
#---------#
TRLPart1 = rfTRL(L=0.435, OD=0.041, ID=0.017, dx=500)
TRLPart2 = rfTRL(L=1.26, OD=0.041, ID=0.017, dx=500)
TRLPart3 = rfTRL(L=1.505, OD=0.041, ID=0.017, dx=500)
TRLPart4 = rfTRL(L=1.26, OD=0.041, ID=0.017, dx=500)
TRLPart5 = rfTRL(L=0.435, OD=0.041, ID=0.017, dx=500)
#load S matrix of antenna
Antenna = rfObject(touchstone='tomas_icrh_linear_2017-vacuum.s2p')
RLCLeft = rfRLC(Cs=CsVal,Cp=CpVals[0])
CaRight = rfRLC(Cs=CaVal)

ct = rfCircuit(Id="MatchingCircuit")

ct.addblock('TL1',TRLPart1,ports=['A','B'],relpos=0)
ct.addblock('TL2',TRLPart2,ports=['B','C'],relpos=TRLPart1.L)
ct.addblock('TL3',TRLPart3,ports=['C','D'],relpos=TRLPart1.L + TRLPart2.L)
ct.addblock('TL4',TRLPart4,ports=['D','E'],relpos=TRLPart1.L + TRLPart2.L + TRLPart3.L)
ct.addblock('TL5',TRLPart5,ports=['E','F'],relpos=TRLPart1.L + TRLPart2.L + TRLPart3.L + TRLPart4.L)

EndOfLine =TRLPart1.L + TRLPart2.L + TRLPart3.L + TRLPart4.L + TRLPart5.L

ct.addblock('Matching',RLCLeft,ports=['F','G'])
ct.addblock('RLoss',rfRLC(Rs=1),ports=['in','out'])
ct.addblock('Antenna',Antenna,ports=['G','H'])
ct.addblock('PreMatch',CaRight,ports=['H','I'])

ct.connect('TL1.A','Source')
ct.connect('TL1.B','TL2.B')
ct.connect('TL2.C','TL3.C')
ct.connect('TL3.D','TL4.D')
ct.connect('TL4.E','TL5.E')
ct.connect('TL5.F','Matching.F')
ct.connect('Matching.G','RLoss.in')
ct.connect('RLoss.out','Antenna.G')
ct.connect('Antenna.H','PreMatch.H')
ct.terminate('PreMatch.I',Z=0) #Grounding

for CpVal in CpVals:
    ct.set('Matching.Cp',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*6*1E3)}, nodes=['TL1.A','TL2.B','TL3.C','TL4.D','TL5.F','PreMatch.D'])

    """
    Calculate error values
    """
    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([0.235,0.66,0.66+0.235+0.795,0.235+0.66+0.795+0.66])
    voltages = np.zeros(len(MeasurePoints))
    Vs = np.zeros(len(MeasurePoints))

    #constants:
    beta = 2*np.pi*FREQ/(3*(10**8))
    S = np.sin(2*beta*MeasurePoints) #array
    C = np.cos(2*beta*MeasurePoints) #array
    BigS = (S[0] - S[1])/(S[2] - S[3])
    BigC = (C[0] - C[1])/(C[2] - C[3])

    #Make array of voltages on the various points:
    voltages[0] = np.abs(Solution['TL2.B'][0])
    voltages[1] = np.abs(Solution['TL3.C'][0])
    voltages[2] = np.abs(Solution['TL4.D'][0])
    voltages[3] = np.abs(Solution['TL5.E'][0])

    Vf = Solution['TL1.A'][2]
    Vr = Solution['TL1.A'][3]

    Vs = (voltages**2)/(np.abs(Vf)**2) - 1
    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC)
    EpsB = -2*v/((1+u)**2 + v**2)
    EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))

    PlotEpsB.append(EpsB)
    PlotEpsG.append(EpsG)

    
plt.plot(CpVals,PlotEpsB,label="$\epsilon_B$")
plt.plot(CpVals,PlotEpsG,label="$\epsilon_G$")
plt.ylabel("$\epsilon$")
plt.xlabel("step")
plt.legend()
plt.show()
