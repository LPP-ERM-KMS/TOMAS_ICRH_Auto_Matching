import matplotlib.pyplot as plt
import numpy as np


thetas = np.linspace(0,2*np.pi,1000)
r = 1

RGammas = r*np.cos(thetas)
CGammas = r*np.sin(thetas)

Gammas = RGammas + CGammas*1j

Ys = (1-Gammas)/(1+Gammas)

#plt.plot(Gammas.real,Gammas.imag)
plt.plot(Ys.real,Ys.imag)
plt.show()
