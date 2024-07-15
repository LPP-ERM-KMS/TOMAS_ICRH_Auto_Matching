import numpy as np
import matplotlib.pyplot as plt

def wave(Gamma,beta):
    return np.abs(np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta))

fs = np.linspace(25e6,55e6,100)
c = 3e8

x = np.linspace(0,2.59,1000)
VoltMeterPositions = np.array([0.235,0.895,1.69,2.35]) 

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

        BigS = (S[0] - S[1])/(S[2] - S[3])
        BigC = (C[0] - C[1])/(C[2] - C[3])

        Vs = (V**2)/(np.abs(Vf)**2) - 1
        print(Vs)

        u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
        v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC)
        EpsB = 2*v/((1+u)**2 + v**2)
        EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))


        plt.quiver(i,j,EpsB,EpsG)

plt.title(f"freq = {f/1e6}MHz")
plt.show()
