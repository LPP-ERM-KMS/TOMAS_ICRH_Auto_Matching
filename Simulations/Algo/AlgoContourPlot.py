import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
four4 = True

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
#x = np.linspace(0,6,1000)
#VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0])*6/2.59
X = []
W = []
Z = []

contourres = 200
f = 25e6
fun = np.zeros((contourres,contourres))

for m,i in enumerate(np.linspace(-1,1,contourres)):
    for n,j in enumerate(np.linspace(-1,1,contourres)):
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

        Vs = (np.abs(V)**2)/(np.abs(Vf)**2) - 1

        S = np.sin(2*beta*VoltMeterPositions) #array
        C = np.cos(2*beta*VoltMeterPositions) #array

        BigS = (S[0] - S[1])/(S[2] - S[3])
        BigC = (C[0] - C[1])/(C[2] - C[3])

        if four4:
            u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[2] - C[3])*(BigC - BigS))
            v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[2] - S[3])*(BigC - BigS))
        else:
            u = (1/2)*(Vs[0]*(1 - S[1]) - Vs[1]*(1 - S[0]))/( (1+C[0])*(1-S[1]) - (1+C[1])*(1-S[0]) )
            v = (-1/2)*(Vs[0]*(1 - C[1]) - Vs[1]*(1 - C[0]))/( (1+C[0])*(1-S[1]) - (1+C[1])*(1-S[0]) )

        RYCorr = 1 - (1 - u**2 - v**2)/((1+u)**2 + v**2)
        IYCorr = 2*v/((1+u)**2 + v**2)

        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsG = NewGamma.real - Gamma.real 
        EpsA = NewGamma.imag - Gamma.imag
        X.append(Gamma.real)
        W.append(Gamma.imag)
        Z.append(NewGamma)
        fun[m,n] = EpsG**2 + EpsA**2

plt.contour(np.linspace(-1,1,contourres),np.linspace(-1,1,contourres),fun)
plt.show()
