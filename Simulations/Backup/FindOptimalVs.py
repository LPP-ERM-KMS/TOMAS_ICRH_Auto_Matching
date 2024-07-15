import numpy as np
import matplotlib.pyplot as plt

def wave(Gamma,beta):
    return np.abs(np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta))

fs = np.linspace(25e6,55e6,100)
c = 3e8


x = np.linspace(0,2.59,1000)
VoltMeterPositions = np.array([0.235,0.895,1.69,2.35])
Gamma = 0.3 + 0.0*1j
V0V1 = np.zeros(len(fs))
V1V2 = np.zeros(len(fs))
V2V3 = np.zeros(len(fs))
V1V3 = np.zeros(len(fs))
V0V2 = np.zeros(len(fs))
V0V3 = np.zeros(len(fs))
EpsAV0V3 = np.zeros(len(fs))
EpsAV0V2 = np.zeros(len(fs))
EpsGV0V3 = np.zeros(len(fs))
EpsGV0V2 = np.zeros(len(fs))

EpsAV1V3 = np.zeros(len(fs))
EpsAV1V2 = np.zeros(len(fs))
EpsGV1V3 = np.zeros(len(fs))
EpsGV1V2 = np.zeros(len(fs))


for j,f in enumerate(fs):
    lam = f/c
    beta = 2*np.pi*lam
    V = np.zeros(4)
    for i,VoltMeterPosition in enumerate(VoltMeterPositions):
        V[i] = wave(Gamma,beta)[np.abs(x-VoltMeterPosition).argmin()]

    Vf = 1
    Vr = Gamma

    S = np.sin(2*beta*VoltMeterPositions) #array
    C = np.cos(2*beta*VoltMeterPositions) #array
    

    EpsAV0V3[j] = (V[0] - Vf)*S[3] - (V[3] - Vf)*S[0]
    EpsGV0V3[j] = (V[0] - Vf)*C[3] - (V[3] - Vf)*C[0]

    EpsAV0V2[j] = (V[0] - Vf)*S[2] - (V[2] - Vf)*S[0]
    EpsGV0V2[j] = (V[0] - Vf)*C[2] - (V[2] - Vf)*C[0]

    EpsAV1V3[j] = (V[1] - Vf)*S[3] - (V[3] - Vf)*S[1]
    EpsGV1V3[j] = (V[1] - Vf)*C[3] - (V[3] - Vf)*C[1]

    EpsAV1V2[j] = (V[1] - Vf)*S[2] - (V[2] - Vf)*S[1]
    EpsGV1V2[j] = (V[1] - Vf)*C[2] - (V[2] - Vf)*C[1]

    V0V1[j] = abs(V[0] - V[1])
    V1V2[j] = abs(V[1] - V[2])
    V2V3[j] = abs(V[2] - V[3])
    V1V3[j] = abs(V[1] - V[3])
    V0V2[j] = abs(V[0] - V[2])
    V0V3[j] = abs(V[0] - V[3])
    
plt.plot(fs/1e6,V0V1,label="V0-V1")
plt.plot(fs/1e6,V1V2,label="V1-V2")
plt.plot(fs/1e6,V2V3,label="V2-V3")
plt.plot(fs/1e6,V1V3,label="V1-V3")
plt.plot(fs/1e6,V0V2,label="V0-V2")
plt.plot(fs/1e6,V0V3,label="V0-V3")
plt.xlabel("MHz")
plt.ylabel("Voltage")
plt.legend()
plt.show()
"""
plt.plot(fs/1e6,EpsAV0V3,label="V0 and V3")
plt.plot(fs/1e6,EpsGV0V3,label="V0 and V3")
plt.plot(fs/1e6,EpsAV0V2,label="V0 and V2")
plt.plot(fs/1e6,EpsGV0V2,label="V0 and V2")
plt.xlabel("MHz")
plt.ylabel("EpsilonG and EpsilonA")
plt.legend()
plt.show()
"""
plt.plot(fs/1e6,np.abs(EpsAV0V3) + np.abs(EpsGV0V3),label="V0 and V3")
plt.plot(fs/1e6,np.abs(EpsAV0V2) + np.abs(EpsGV0V2),label="V0 and V2")
plt.plot(fs/1e6,np.abs(EpsAV1V3) + np.abs(EpsGV1V3),label="V1 and V3")
plt.plot(fs/1e6,np.abs(EpsAV1V2) + np.abs(EpsGV1V2),label="V1 and V2")
plt.xlabel("MHz")
plt.ylabel("Epsilon")
plt.legend()
plt.show()


plt.plot(fs/1e6,EpsAV1V2,label="V1 and V2")
plt.plot(fs/1e6,EpsGV1V2,label="V1 and V2")
plt.xlabel("MHz")
plt.ylabel("Epsilon")
plt.legend()
plt.show()
