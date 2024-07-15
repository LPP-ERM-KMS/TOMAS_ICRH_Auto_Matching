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

fig, axs = plt.subplots(2,3,sharex=True,sharey=True)

x = np.linspace(0,2.59,1000)
VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0]) #added V+/V- as bonus measurement 
#x = np.linspace(0,6,1000)
#VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0])*6/2.59
V1V2Errors = []
V1V3Errors = []
V0V3Errors = []
V0V1Errors = []
V0V2Errors = []
V2V3Errors = []

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
        

        IYCorr = (V[0] - Vf)*S[3] - (V[3] - Vf)*S[0]
        RYCorr = (V[0] - Vf)*C[3] - (V[3] - Vf)*C[0]

        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsGV0V3 = NewGamma.real - Gamma.real
        EpsAV0V3 = NewGamma.imag - Gamma.imag

        V0V3Errors.append([i,j,EpsGV0V3,EpsAV0V3])

        IYCorr = (V[0] - Vf)*S[1] - (V[1] - Vf)*S[0]
        RYCorr = (V[0] - Vf)*C[1] - (V[1] - Vf)*C[0]

        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsGV0V1 = NewGamma.real - Gamma.real
        EpsAV0V1 = NewGamma.imag - Gamma.imag
        
        V0V1Errors.append([i,j,EpsGV0V1,EpsAV0V1])


        IYCorr = (V[0] - Vf)*S[2] - (V[2] - Vf)*S[0]
        RYCorr = (V[0] - Vf)*C[2] - (V[2] - Vf)*C[0]

        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsGV0V2 = NewGamma.real - Gamma.real
        EpsAV0V2 = NewGamma.imag - Gamma.imag

        V0V2Errors.append([i,j,EpsGV0V2,EpsAV0V2])
       
        IYCorr = (V[1] - Vf)*S[3] - (V[3] - Vf)*S[1]
        RYCorr = (V[1] - Vf)*C[3] - (V[3] - Vf)*C[1]

        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsGV1V3 = NewGamma.real - Gamma.real
        EpsAV1V3 = NewGamma.imag - Gamma.imag

        V1V3Errors.append([i,j,EpsGV1V3,EpsAV1V3])

        IYCorr = (V[1] - Vf)*S[2] - (V[2] - Vf)*S[1]
        RYCorr = (V[1] - Vf)*C[2] - (V[2] - Vf)*C[1]

        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsGV1V2 = NewGamma.real - Gamma.real
        EpsAV1V2 = NewGamma.imag - Gamma.imag

        V1V2Errors.append([i,j,EpsGV1V2,EpsAV1V2])

        IYCorr = (V[2] - Vf)*S[3] - (V[3] - Vf)*S[2]
        RYCorr = (V[2] - Vf)*C[3] - (V[3] - Vf)*C[2]
        
        Y_n = Y + RYCorr + 1j*IYCorr
        NewGamma = (1 - Y_n)/(1 + Y_n)
        EpsGV2V3 = NewGamma.real - Gamma.real
        EpsAV2V3 = NewGamma.imag - Gamma.imag

        V2V3Errors.append([i,j,EpsGV2V3,EpsAV2V3])

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

InitialPoint = 0.4+1j*0.3

V1V2Errors = np.array(V1V2Errors)
V1V3Errors = np.array(V1V3Errors)
V0V3Errors = np.array(V0V3Errors)
V0V1Errors = np.array(V0V1Errors)
V0V2Errors = np.array(V0V2Errors)
V2V3Errors = np.array(V2V3Errors)

Route = FindRoute(InitialPoint,V0V1Errors)
for point in Route:
    axs[0,0].scatter(point[0],point[1],color="blue")

Route = FindRoute(InitialPoint,V0V2Errors)
for point in Route:
    axs[1,0].scatter(point[0],point[1],color="blue")

Route = FindRoute(InitialPoint,V1V3Errors)
for point in Route:
    axs[1,1].scatter(point[0],point[1],color="blue")

Route = FindRoute(InitialPoint,V1V2Errors)
for point in Route:
    axs[0,1].scatter(point[0],point[1],color="blue")

Route = FindRoute(InitialPoint,V0V3Errors)
for point in Route:
    axs[1,2].scatter(point[0],point[1],color="blue")

Route = FindRoute(InitialPoint,V2V3Errors)
for point in Route:
    axs[0,2].scatter(point[0],point[1],color="blue")

#plt.title(f"freq = {f/1e6}MHz")
fig.text(0.5, 0.04, '$\Re(\Gamma)$', ha='center')
fig.text(0.04, 0.5, '$\Im(\Gamma)$', va='center', rotation='vertical')
fig.suptitle(f'Frequency = {f}MHz', fontsize=16)
plt.show()


