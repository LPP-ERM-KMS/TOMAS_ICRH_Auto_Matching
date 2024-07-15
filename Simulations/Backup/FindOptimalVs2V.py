import numpy as np
import matplotlib.pyplot as plt

def wave(Gamma,beta):
    return np.abs(np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta))

fs = np.linspace(25e6,55e6,100)
c = 3e8

fig, axs = plt.subplots(2,3)

x = np.linspace(0,2.59,1000)
VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0]) #added V+/V- as bonus measurement 
#x = np.linspace(0,6,1000)
#VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0])*6/2.59

f = 25e6

for i in np.linspace(-1,1,20):
    for j in np.linspace(-1,1,20):
        Gamma = i + 1j*j
        lam = f/c
        beta = 2*np.pi*lam
        V = np.zeros(len(VoltMeterPositions))
        for k,VoltMeterPosition in enumerate(VoltMeterPositions):
            V[k] = wave(Gamma,beta)[np.abs(x-VoltMeterPosition).argmin()]

        Vf = 1

        S = np.sin(2*beta*VoltMeterPositions) #array
        C = np.cos(2*beta*VoltMeterPositions) #array
        

        EpsAV0V3 = (V[0] - Vf)*S[3] - (V[3] - Vf)*S[0]
        EpsGV0V3 = (V[0] - Vf)*C[3] - (V[3] - Vf)*C[0]

        EpsAV0V1 = (V[0] - Vf)*S[1] - (V[1] - Vf)*S[0]
        EpsGV0V1 = (V[0] - Vf)*C[1] - (V[1] - Vf)*C[0]

        EpsAV0V2 = (V[0] - Vf)*S[2] - (V[2] - Vf)*S[0]
        EpsGV0V2 = (V[0] - Vf)*C[2] - (V[2] - Vf)*C[0]

        EpsAV1V3 = (V[1] - Vf)*S[3] - (V[3] - Vf)*S[1]
        EpsGV1V3 = (V[1] - Vf)*C[3] - (V[3] - Vf)*C[1]

        EpsAV1V2 = (V[1] - Vf)*S[2] - (V[2] - Vf)*S[1]
        EpsGV1V2 = (V[1] - Vf)*C[2] - (V[2] - Vf)*C[1]

        EpsAV2V3 = (V[2] - Vf)*S[3] - (V[3] - Vf)*S[2]
        EpsGV2V3 = (V[2] - Vf)*C[3] - (V[3] - Vf)*C[2]

        EpsAV0V4 = (V[0] - Vf)*S[4] - (V[4] - Vf)*S[0]
        EpsGV0V4 = (V[0] - Vf)*C[4] - (V[4] - Vf)*C[0]

        EpsAV2V4 = (V[2] - Vf)*S[4] - (V[4] - Vf)*S[2]
        EpsGV2V4 = (V[2] - Vf)*C[4] - (V[4] - Vf)*C[2]

        EpsAV3V4 = (V[3] - Vf)*S[4] - (V[4] - Vf)*S[3]
        EpsGV3V4 = (V[3] - Vf)*C[4] - (V[4] - Vf)*C[3]


        #plt.quiver(i,j,EpsAV1V2,EpsGV1V2)
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

#        axs[0,2].quiver(i,j,EpsGV2V4,EpsAV2V4)
#        axs[1,2].quiver(i,j,EpsGV3V4,EpsAV3V4)
#        axs[2,2].quiver(i,j,EpsGV0V4,EpsAV0V4)

plt.title(f"freq = {f/1e6}MHz")
plt.show()
