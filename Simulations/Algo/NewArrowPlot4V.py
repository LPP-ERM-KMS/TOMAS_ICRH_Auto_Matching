import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt

def wave(Gamma,beta):
    return np.abs(np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta))

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
#x = np.linspace(0,6,1000)
#VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0])*6/2.59

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

        S0123 = (S[0] - S[1])/(S[2] - S[3])
        C0123 = (C[0] - C[1])/(C[2] - C[3])
        S2301 = (S[2] - S[3])/(S[0] - S[1])
        C2301 = (C[2] - C[3])/(C[0] - C[1])

        Vs = (np.abs(V)**2/np.abs(Vf)**2) - 1

        udenom = 2*(C[2] - C[3])*(C0123 - S0123)
        vdenom = 2*(S[0] - S[1])*(C2301 - S2301)

        u = ((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*S0123)/udenom
        v = ((Vs[2] - Vs[3]) - (Vs[0] - Vs[1])*C2301)/vdenom

        RYCorr = 2*v/((1+u)**2 + v**2)
        IYCorr = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))

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

for i in np.linspace(-2,2,30):
    for j in np.linspace(-2,2,30):
        Y = i + 1j*j
        Gamma = (1-Y)/(1+Y)
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

        Vs = (np.abs(V)**2)/(np.abs(Vf)**2) 

        u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
        v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

        IYCorr = 2*v/((1+u)**2 + v**2)
        RYCorr = (1 - ((1-u**2-v**2)/((1+u)**2 + v**2)))

        plt.quiver(i,j,RYCorr,IYCorr)

plt.title("New Algorithm")
plt.xlabel("$\Re(y))$")
plt.ylabel("$\Im(y))$")
plt.show()

