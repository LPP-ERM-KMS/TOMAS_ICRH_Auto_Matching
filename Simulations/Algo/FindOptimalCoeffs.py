import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod as TomasCircuit
from alive_progress import alive_bar

MHz_ = 1e6
pF_ = 1e-12

CsVals = np.linspace(120e-12,170e-12,50)#capacitor is defined in steps
CpVals = np.linspace(100e-12,750e-12,200)#capacitor is defined in steps

U = []
V = []


def Algo4V(Solution):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([0.235,0.895,1.69,2.35])
    V = np.zeros(len(MeasurePoints))

    #constants:
    beta = 2*np.pi*FREQ/(3*(10**8))
    S = np.sin(2*beta*MeasurePoints) #array
    C = np.cos(2*beta*MeasurePoints) #array
    BigS = (S[0] - S[1])/(S[2] - S[3])
    BigC = (C[0] - C[1])/(C[2] - C[3])


    #Make array of voltages on the various points:
    V[3] = np.abs(Solution['V1toV0.V0'][0])
    V[2] = np.abs(Solution['V1toV0.V1'][0])
    V[1] = np.abs(Solution['V3toV2.V2'][0])
    V[0] = np.abs(Solution['V3toV2.V3'][0]) 

    Vf = abs(Solution['V3toV2.V3'][2])
    Vr = abs(Solution['V3toV2.V3'][3])

    Gamma = np.abs(Vr/Vf)


    Vs = (V/Vf)**2

    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

    EpsB = 2*v/((1+u)**2 + v**2)
    EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))

    return EpsG,EpsB

def Algo2V(Solution):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([0.235,0.895,1.69,2.35])
    V = np.zeros(len(MeasurePoints))

    #constants:
    beta = 2*np.pi*FREQ/(3*(10**8))
    S = np.sin(2*beta*MeasurePoints) #array
    C = np.cos(2*beta*MeasurePoints) #array
    BigS = (S[0] - S[1])/(S[2] - S[3])
    BigC = (C[0] - C[1])/(C[2] - C[3])


    #Make array of voltages on the various points:
    V[0] = np.abs(Solution['V1toV0.V0'][0])
    V[1] = np.abs(Solution['V1toV0.V1'][0])
    V[2] = np.abs(Solution['V3toV2.V2'][0])
    V[3] = np.abs(Solution['V3toV2.V3'][0])

    Vf = abs(Solution['V3toV2.V3'][2])
    Vr = abs(Solution['V3toV2.V3'][3])

    Gamma = np.abs(Vr/Vf)

    EpsB = (V[0] - Vf)*S[3] - (V[3] - Vf)*S[0]
    EpsG = (V[0] - Vf)*C[3] - (V[3] - Vf)*C[0]

    return -1*EpsG,EpsB

FREQ = 25*MHz_

Gammas = []

ct = TomasCircuit()

CaVal = 220*pF_
ct.set('Ca.Cs',CaVal)

with alive_bar(len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)

            Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V3toV2'])

            Vf = abs(Solution['V3toV2.V3'][2])
            Vr = abs(Solution['V3toV2.V3'][3])

            #epsg,epsb = Algo4V(Solution)
            epsg,epsb = Algo2V(Solution)
            x = np.array([epsg,epsb])
            x = 2*x / np.linalg.norm(x)
            epsg = x[0]
            epsb = x[1]

            U.append(epsb)
            V.append(epsg)

            
            Gammas.append(Vr/Vf)
            bar()

    X, Y = np.meshgrid(CpVals/pF_,CsVals/pF_) 

    #X has form [[CpVals],[CpVals],..]
    #Y has form [[CsVals[0]],[CsVals[1]],..]
    #I.e a MESHGRID is created from the x and y axes (as the name suggests)

    Z = np.array(Gammas).reshape(len(CsVals),len(CpVals)) #put into shape with x length CsVals and y length CpVals
    levels = np.linspace(0.0, 1.0, 11)
    CS = plt.contourf(X, Y, Z, levels=levels, cmap=cm.coolwarm, extend='min')
    colorbar = plt.colorbar(CS)
    plt.ylabel("Cs (pF)")
    plt.xlabel("Cp (pF)")
    plt.title(f"Ca = {round(CaVal/pF_,2)}, Freq= {FREQ/MHz_}MHz; quiver shows the 2V textor algorithm")

U = np.array(U)
V = np.array(V)
# Normalize the arrows:
#U = U / np.sqrt(U**2 + V**2);
#V = V / np.sqrt(U**2 + V**2);

plt.quiver(X,Y,U,V)

plt.show()
