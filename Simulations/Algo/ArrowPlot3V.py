import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt

def wave(Gamma,beta):
    return np.abs(np.exp(-1j*x*beta) + Gamma*np.exp(1j*x*beta))

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

x = np.linspace(0,2.59,1000)
VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0]) #added V+/V- as bonus measurement 

f = 25e6

for i in np.linspace(-1,1,30):
    for j in np.linspace(-1,1,30):
        if i**2 + j**2 > 1:
            continue
        Gamma = i + 1j*j
        Y = (1-Gamma)/(1+Gamma)
        lam = f/c
        beta = 2*np.pi*lam
        V = np.zeros(len(VoltMeterPositions))
        for k,VoltMeterPosition in enumerate(VoltMeterPositions):
            V[k] = wave(Gamma,beta)[np.abs(x-VoltMeterPosition).argmin()]

        Vf = 1

        S = np.sin(2*beta*VoltMeterPositions) #array
        C = np.cos(2*beta*VoltMeterPositions) #array

        Ai = 0
        Bi = 2
        Ci = 3

        BigS = (S[0] - S[1])/(S[1] - S[2])
        BigC = (C[0] - C[1])/(C[1] - C[2])
        SABBC = (S[Ai] - S[Bi])/(S[Bi] - S[Ci])
        CABBC = (C[Ai] - C[Bi])/(C[Bi] - C[Ci])

       
        Vs = (V/Vf)**2 - 1

        udenom = 2*(C[Bi]-C[Ci])*(CABBC - SABBC)
        vdenom = 2*(S[Bi]-S[Ci])*(CABBC - SABBC)

        u = ((Vs[Ai] - Vs[Bi]) - (Vs[Bi] - Vs[Ci])*SABBC)/udenom
        v = ((Vs[Ai] - Vs[Bi]) - (Vs[Bi] - Vs[Ci])*CABBC)/vdenom
        
        IYCorr = -2*v/((1+u)**2 + v**2)
        RYCorr = ((1-u**2-v**2)/((1+u)**2 + v**2)) - 1

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


