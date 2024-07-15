import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs

def findClosestValue(givenList, target):
    def difference(givenList):
        return abs(givenList - target)
    result = min(givenList, key=difference)
    return result


CsVals = np.linspace(25e-12,1000e-12,1000)#capacitor is defined in steps
CpVals = np.linspace(7e-12,1000e-12,1000)#capacitor is defined in steps

MHz_ = 1e6
pF_ = 1e-12

FREQ = 25*MHz_

CaVal=150*pF_ #CalcCa(FREQ*(10**(-6)))*(10**(-12))

#initial values:
CsVal = 210*pF_ #CsVals[100]
CpVal = 380*pF_ #CpVals[100]
CsFactor = 0.000001*pF_
CpFactor = 0.000001*pF_

PlotCs = []
PlotCp = []
admittances = []
ReflCoeffs = []

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
Zbase = 50
ct = rfCircuit(Zbase = Zbase)

# antenna strap

A2Ca1_L = 0.1715
A2Ca1_Z = 50
T2Cs1_L = 0.121
T2Cs1_Z = 50

tsStrap = 'Antenna/Tomas-Ref_geo-R=200-Diel_eps=0500.s2p'


strap = rfObject(touchstone=tsStrap, ports=['cap','t'])

if strap.fs[0] == 0.:
    strap.fs[0] = 1e-3 # avoid divisions by 0 at 0 Hz in deembed

strap.deembed({'cap':(A2Ca1_L, A2Ca1_Z), 
                    't'  :(T2Cs1_L, T2Cs1_Z)})
        
# circuit
T2Cs_L = 0.16
T2Cs_Z = 50
Cs2Cp_L = 0.2
Cs2Cp_Z = 50
A2Ca_L = 0.29
A2Ca_Z = 50.
CptoV0_Z = 50
V0toV1_Z = 50
V1toV2_Z = 50
V2toV3_Z = 50

CptoV0_L = 0.235
V0toV1_L = 0.66
V1toV2_L = 0.795
V2toV3_L = 0.66

LsCaps_H = 30e-9

# build the circuit
ct.addblock('strap', strap, 
                 # ports=['t','cap']
                 )
ct.addblock('A2Ca1',
    rfTRL(L=A2Ca1_L, Z0TL=A2Ca1_Z,
    ports= ['ant','vcw'])
)
ct.connect('strap.cap','A2Ca1.ant')
ct.addblock('A2Ca', 
    rfTRL(L=A2Ca_L, Z0TL=A2Ca_Z,
    ports= ['vcw','cap'])
)
ct.connect('A2Ca1.vcw', 'A2Ca.vcw')
ct.addblock('Ca', rfRLC(Cs=CaVal, Ls=LsCaps_H, ports= ['cap','sc']))
#ct.addblock('RLoss', rfRLC(Rs=0.0, ports= ['in','out']))
ct.connect('Ca.cap', 'A2Ca.cap')
#ct.connect('RLoss.in', 'A2Ca.cap')
ct.terminate('Ca.sc',Z=0.)
ct.addblock('T2Cs1',
    rfTRL(L=T2Cs1_L, Z0TL=T2Cs1_Z,
    ports= ['t','vcw'])
)
ct.connect('strap.t', 'T2Cs1.t')
ct.addblock('T2Cs',
    rfTRL(L=T2Cs_L, Z0TL=T2Cs_Z,
    ports= ['vcw','cap'])
)
ct.connect('T2Cs1.vcw', 'T2Cs.vcw')
ct.addblock('Cs', rfRLC(Cs=CsVal, Ls=LsCaps_H, ports= ['cap','toCp']))
ct.connect('T2Cs.cap','Cs.cap')
ct.addblock('Cs2Cp',
    rfTRL(L=Cs2Cp_L, Z0TL=Cs2Cp_Z,
    ports= ['Cs','Cp'])
)
ct.connect('Cs.toCp', 'Cs2Cp.Cs')
ct.addblock('Cp', rfRLC(Cs=CpVal, Ls=LsCaps_H, ports= ['Cp','sc']))
ct.terminate('Cp.sc', Z=0)


# Coaxial cable
ct.addblock('CptoV0', 
    rfTRL(L=CptoV0_L, Z0TL=CptoV0_Z,
    ports= ['Cp','V0'])
)

ct.connect('Cs2Cp.Cp', 'Cp.Cp', 'CptoV0.Cp')

ct.addblock('V0toV1', 
    rfTRL(L=V0toV1_L, Z0TL=V0toV1_Z,
    ports= ['V0','V1'])
)
ct.addblock('V1toV2', 
    rfTRL(L=V1toV2_L, Z0TL=V1toV2_Z,
    ports= ['V1','V2'])
)
ct.connect('CptoV0.V0','V0toV1.V0')
ct.connect('V0toV1.V1','V1toV2.V1')

ct.addblock('V2toV3', 
    rfTRL(L=V2toV3_L, Z0TL=V2toV3_Z,
    ports= ['V2','V3'])
)
ct.connect('V1toV2.V2','V2toV3.V2')
ct.connect('V2toV3.V3','Source')

while not Converged:
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

#    Calculate new values of Cs and Cp

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
    ReflCoeffs.append(ReflCoeff)

    EpsB = (V[1] - Vf)*S[3] - (V[3] - Vf)*S[1]
    EpsG = (V[1] - Vf)*C[3] - (V[3] - Vf)*C[1]

    print("Eps:")
    print(EpsB + 1j*EpsG)
    admittance = EpsG+1 + 1j*EpsB
    admittances.append(admittance)
    
    DeltaCs = CsFactor*EpsG
    DeltaCp = CpFactor*EpsB
    #gm = ct.getS(25e6)[0]
    #tSol = ct.Solution(25e6, {'Source':np.sqrt(2*50*5*1E3)}, nodes=['CptoV0.Cp'])
    #yCp = -tSol['CptoV0.Cp'][1]/tSol['CptoV0.Cp'][0]*50
    #DeltaCs = max(-10E-12, min((yCp.real-1), 10E-12))
    #DeltaCp = max(-10E-12, min(yCp.imag*(-1), 10E-12))


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
        CpFactor *= -0.8
    if step > 500:
        maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': np.sqrt(2*50*6*1E3)})
        break


print("Final capacitor values:")
print(f"Cs={CsVal} and Cp={CpVal}")
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
