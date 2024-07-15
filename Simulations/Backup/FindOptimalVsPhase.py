import numpy as np
import matplotlib.pyplot as plt

def wave(Gamma,beta):
    return np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta)
def UandV(VBar,Vtilde,C,S,Probe1,Probe2):
    u = (C[Probe2]*VBar[Probe1] - S[Probe1]*Vtilde[Probe2])/(C[Probe1]*C[Probe2] + S[Probe1]*S[Probe1])
    v = (S[Probe1]*VBar[Probe1] + C[Probe1]*Vtilde[Probe2])/(C[Probe1]*C[Probe2] + S[Probe1]*S[Probe1])
    return u,v
"""
def UandV(VBar,Vtilde,C,S,Probe1,Probe2):
    u = (VBar[Probe1] - S[Probe1]*VBar[Probe2])/(C[Probe1]*S[Probe2] - S[Probe1]*C[Probe2])
    v = (Vtilde[Probe1] - S[Probe1]*Vtilde[Probe2])/(C[Probe1]*S[Probe2] - S[Probe1]*C[Probe2])
    return u,v
"""
def Errors(u,v):
    EpsR = 1 - (1 - u**2 - v**2)/((1+u)**2 + v**2)
    EpsC = -2*v/((1+u)**2 + v**2)
    return EpsR,EpsC

fs = np.linspace(25e6,55e6,100)
c = 3e8

#fig, axs = plt.subplots(2,3)

x = np.linspace(0,2.59,1000)
VoltMeterPositions = np.array([0.235,0.895,1.69,2.35])
#x = np.linspace(0,6,1000)
#VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0])*6/2.59

f = 25e6

for i in np.linspace(-1,1,20):
    for j in np.linspace(-1,1,20):
        if i**2 + j**2 > 1:
            continue
        Gamma = i + 1j*j
        y = (1-Gamma)/(1+Gamma)
        lam = f/c
        beta = 2*np.pi*lam
        V = np.zeros(len(VoltMeterPositions))
        VBar = np.zeros(len(VoltMeterPositions))
        Vtilde = np.zeros(len(VoltMeterPositions))
        for k,VoltMeterPosition in enumerate(VoltMeterPositions):
            V[k] = wave(Gamma,beta)[np.abs(x-VoltMeterPosition).argmin()]

        S = np.sin(2*beta*VoltMeterPositions) #array
        C = np.cos(2*beta*VoltMeterPositions) #array


        VBar = V.real - C
        Vtilde = V.imag - S

        #V1 and V3
        u,v = UandV(VBar,Vtilde,C,S,1,2)
        EpsRV1V3, EpsCV1V3 = Errors(u,v)
      
        """
        axs[0,0].quiver(i,j,EpsGV0V1,EpsAV0V1)
        axs[0,0].set_title("Error values using V0 and V1")
        axs[0,1].quiver(i,j,EpsGV1V2,EpsAV1V2)
        axs[0,1].set_title("Error values using V1 and V2")
        axs[0,2].quiver(i,j,EpsGV2V3,EpsAV2V3)
        axs[0,2].set_title("Error values using V2 and V3")

        axs[1,0].quiver(i,j,EpsGV0V2,EpsAV0V2)
        axs[1,0].set_title("Error values using V0 and V2")
        axs[1,1].quiver(i,j,EpsGV1V3,EpsAV1V3)
        axs[1,1].set_title("Error values using V1 and V3")
        axs[1,2].quiver(i,j,EpsGV0V3,EpsAV0V3)
        axs[1,2].set_title("Error values using V0 and V3")
        """
        plt.quiver(i,j,EpsRV1V3,EpsCV1V3)

plt.show()
