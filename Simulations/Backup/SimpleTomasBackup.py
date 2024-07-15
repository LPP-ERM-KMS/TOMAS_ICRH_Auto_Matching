import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs

CsVals = np.linspace(25e-12,1000e-12,1000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,1000)#capacitor is defined in steps

FREQ = 25*10**6

pF_ = 1e-12

CaVal = 194#CalcCa(FREQ*(10**(-6)))*(10**(-12))

#initial values:
CsVal = 133#CsVals[100]
CpVal = 48#CpVals[100]
CsFactor = 0.01*pF_
CpFactor = 0.01*pF_

PlotCs = []
PlotCp = []
admittances = []
ReflCoeffs = []
ReflCoeffs.append(1)

PlotCs.append(CsVal)
PlotCp.append(CpVal)

PlotEpsB = []
PlotEpsG = []

Steps = []
step = 1
Steps.append(step)
Converged = False

#---------#
# Circuit #
#---------#
TRLPart1 = rfTRL(L=0.235, OD=0.041, ID=0.017, dx=500)
TRLPart2 = rfTRL(L=0.66, OD=0.041, ID=0.017, dx=500)
TRLPart3 = rfTRL(L=0.795, OD=0.041, ID=0.017, dx=500)
TRLPart4 = rfTRL(L=0.66, OD=0.041, ID=0.017, dx=500)
TRLPart5 = rfTRL(L=0.235, OD=0.041, ID=0.017, dx=500)
#load S matrix of antenna
Antenna = rfObject(touchstone='tomas_icrh_linear_2017-vacuum.s2p')
RLCLeft = rfRLC(Cs=CsVal,Cp=CpVal)
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
ct.connect('Matching.G','Antenna.G')
ct.connect('Antenna.H','RLoss.in')
ct.connect('RLoss.out','PreMatch.H')
ct.terminate('PreMatch.I',Z=0) #Grounding

while not Converged:
    ct.set('Matching.Cs',CsVal)
    ct.set('Matching.Cp',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*6*1E3)}, nodes=['TL1.A','TL2.B','TL3.C','TL4.D','TL5.F','PreMatch.D'])

    """
    Calculate new values of Cs and Cp
    """
    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([3.65,3.65+0.66,6.65+0.66+0.795,3.65+0.66+0.795+0.66])
    V = np.zeros(len(MeasurePoints))
    Vs = np.zeros(len(MeasurePoints))

    #constants:
    beta = 2*np.pi*FREQ/(3*(10**8))
    S = np.sin(2*beta*MeasurePoints) #array
    C = np.cos(2*beta*MeasurePoints) #array
    BigS = (S[0] - S[1])/(S[2] - S[3])
    BigC = (C[0] - C[1])/(C[2] - C[3])

    #Make array of voltages on the various points:
    V[0] = np.abs(Solution['TL2.B'][0])
    V[1] = np.abs(Solution['TL3.C'][0])
    V[2] = np.abs(Solution['TL4.D'][0])
    V[3] = np.abs(Solution['TL5.E'][0])

    Vf = abs(Solution['TL1.A'][2])
    Vr = abs(Solution['TL1.A'][3])
    print("reflection coeff:")
    print(np.abs(Vr/Vf))
    ReflCoeff = np.abs(Vr/Vf)
    ReflCoeffs.append(ReflCoeff)
    #Vs = (V**2)/(np.abs(Vf)**2) - 1
    #u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
    #v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC)
    #EpsB = -2*v/((1+u)**2 + v**2)
    #EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))
    EpsB = (V[1] - Vf)*S[3] - (V[3] - Vf)*S[1]
    EpsG = (V[1] - Vf)*C[3] - (V[3] - Vf)*C[1]

    print("Eps:")
    print(EpsB + 1j*EpsG)
    admittance = EpsG+1 + 1j*EpsB
    #Gamma = u + 1j*v
    admittances.append(admittance)
    
    DeltaCs = CsFactor*EpsG
    DeltaCp = CpFactor*EpsB

    CsVal += DeltaCs
    CpVal += DeltaCp
        
    #IndexCloseToNewCs = np.where(findClosestValue(CsVals, CsVal) == CsVals)#get index of closest lying computed 
    #IndexCloseToNewCp = np.where(findClosestValue(CpVals, CpVal) == CpVals)#get index of closest lying computed 

    #CsVal = CsVals[IndexCloseToNewCs][0]
    #CpVal = CpVals[IndexCloseToNewCp][0]

    PlotCs.append(CsVal)
    PlotCp.append(CpVal)

    PlotEpsB.append(EpsB)
    PlotEpsG.append(EpsG)

    #Converged enough:
    if (np.abs(ReflCoeff) < 0.1):
        print("nicely converged")
        maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': np.sqrt(2*50*6*1E3)})
        Converged = True
    step += 1
    Steps.append(step)
    if ((PlotCs[-2] == PlotCs[-1]) and (PlotCp[-2] == PlotCp[-1])):
        print("converged due to a lack of changes")
        maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': np.sqrt(2*50*1E3)})
        Converged = True
    if (ReflCoeffs[-2] < ReflCoeffs[-1]) and (ReflCoeff < 0.2):
        print("overshot")
        CsFactor *= -0.8
        CpFactor *= -0.8


plotVSWs(VSWs) 
plt.show()

plt.plot(np.arange(0,len(ReflCoeffs)),np.real(ReflCoeffs),label="real part")
plt.show()

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
