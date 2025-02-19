import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt

def wave(Gamma,beta):
    return np.abs(np.exp(-1j*x*beta) + Gamma*np.exp(1j*x*beta))
def Vf(Gamma,beta):
    return np.exp(-1j*x*beta)
def Vr(Gamma,beta):
    return Gamma*np.exp(1j*x*beta)

def FindRoute(InitialPoint,V1V2Errors,ScaleFactor=0.010):
    SearchFile = V1V2Errors[:,:2] #search Gamma
    tree = spatial.KDTree(SearchFile)
    Point = InitialPoint
    route = []

    for _ in range(1000):
        route.append([Point.real,Point.imag])
        IndexOfClosest = tree.query(np.array([Point.real,Point.imag]))[1] #get index of closest lying computed yield
        RCorr = V1V2Errors[IndexOfClosest][2]
        CCorr = V1V2Errors[IndexOfClosest][3]
        Point += ScaleFactor*(RCorr + 1j*CCorr)

    return np.array(route)


fs = np.linspace(25e6,55e6,1000)
c = 3e8

x = np.linspace(0,2.585,1000)

f = 25e6

for i in np.linspace(-1,1,30):
    for j in np.linspace(-1,1,30):
        if i**2 + j**2 > 1:
            continue
        Gamma = i + 1j*j
        Y = (1-Gamma)/(1+Gamma)

        lam = f/c
        beta = 2*np.pi*lam
        #Measured at V3
        Vf_ = Vf(Gamma,beta)[np.abs(x-2.35).argmin()]
        Vr_ = Vr(Gamma,beta)[np.abs(x-2.35).argmin()]

        GemetenGamma = Vr_/Vf_

        rho = np.abs(Vr_)/np.abs(Vf_)
        phi =  np.angle(Vr_) - np.angle(Vf_)

        Beta = 2*np.pi*f/c
        BetaL3 = Beta*2.35
        phi = phi - 2*BetaL3 #phase transform from dir coupler to load
        Gamma = rho*np.exp(1j*(phi))

        u = Gamma.real
        v = Gamma.imag

        IYCorr = v
        RYCorr = (u**2 + v**2 + u)

        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsG = NewGamma.real - Gamma.real 
        EpsA = NewGamma.imag - Gamma.imag

        plt.quiver(i,j,EpsG,EpsA)

#plt.title("Algorithm using voltage probes 0-3 at 25MHz")
plt.title("New Algorithm")
plt.xlabel("$\Re(\Gamma))$")
plt.ylabel("$\Im(\Gamma))$")
plt.show()
