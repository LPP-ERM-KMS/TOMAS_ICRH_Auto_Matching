import numpy as np
import matplotlib.pyplot as plt

def wave(Gamma):
    return np.abs(np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta))

fs = [25e6,35e6,55e6]
c = 3e8


x = np.linspace(0,2.59,1000)
Gammas = np.linspace(0,0.5,5)

for f in fs:
    for Gamma in Gammas:
        lam = f/c
        beta = 2*np.pi*lam
        plt.plot(x,wave(Gamma),label="Gamma={0},f={1}".format(Gamma,f))
    
plt.axvline(x=0.235,label="V3")
plt.axvline(x=0.895,label="V2")
plt.axvline(x=1.69,label="V1")
plt.axvline(x=2.35,label="V0")
plt.xlabel("x (m)")
plt.ylabel("V")
plt.legend()
plt.show()
