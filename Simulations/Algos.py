import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod as TomasCircuit

def ModAlgo3V(Solution,FREQ,probeindexA,probeindexB,probeindexC):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    Ai = probeindexA
    Bi = probeindexB
    Ci = probeindexC

    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([2.35,1.69,0.895,0.235])
    V = np.zeros(len(MeasurePoints))

    #constants:
    beta = 2*np.pi*FREQ/(3*(10**8))
    S = np.sin(2*beta*MeasurePoints) #array
    C = np.cos(2*beta*MeasurePoints) #array
    SABBC = (S[Ai] - S[Bi])/(S[Bi] - S[Ci])
    CABBC = (C[Ai] - C[Bi])/(C[Bi] - C[Ci])

    #Make array of voltages on the various points:
    V[0] = np.abs(Solution['V1toV0.V0'][0])
    V[1] = np.abs(Solution['V1toV0.V1'][0])
    V[2] = np.abs(Solution['V3toV2.V2'][0])
    V[3] = np.abs(Solution['V3toV2.V3'][0]) 

    Vf = Solution['V3toV2.V3'][2]
    Vr = Solution['V3toV2.V3'][3]

    Gamma = Vr/Vf

    Vf = np.abs(Vf)
    Vs = (V**2)/(Vf**2) - 1

    udenom = 2*(C[Bi]-C[Ci])*(CABBC - SABBC)
    vdenom = 2*(S[Bi]-S[Ci])*(SABBC - CABBC)

    u = ((Vs[Ai] - Vs[Bi]) - (Vs[Bi] - Vs[Ci])*SABBC)/udenom
    v = ((Vs[Ai] - Vs[Bi]) - (Vs[Bi] - Vs[Ci])*CABBC)/vdenom

    EpsB = 2*v/((1+u)**2 + v**2)
    EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))

    #beta = np.angle(Gamma) #yolo
    beta = -1*np.angle(Gamma) #yolo
    #x = np.array([EpsG-*EpsB,EpsB])
    x = np.array([np.cos(beta)*EpsG-np.sin(beta)*EpsB,np.sin(beta)*EpsG + np.cos(beta)*EpsB])
    x = x/np.linalg.norm(x)
    EpsG = x[0]
    EpsB = x[1]

    return EpsG,EpsB


def Algo3V(Solution,FREQ,probeindexA,probeindexB,probeindexC):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    Ai = probeindexA
    Bi = probeindexB
    Ci = probeindexC

    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([2.35,1.69,0.895,0.235])
    V = np.zeros(len(MeasurePoints))

    #constants:
    beta = 2*np.pi*FREQ/(3*(10**8))
    S = np.sin(2*beta*MeasurePoints) #array
    C = np.cos(2*beta*MeasurePoints) #array
    SABBC = (S[Ai] - S[Bi])/(S[Bi] - S[Ci])
    CABBC = (C[Ai] - C[Bi])/(C[Bi] - C[Ci])

    #Make array of voltages on the various points:
    V[0] = np.abs(Solution['V1toV0.V0'][0])
    V[1] = np.abs(Solution['V1toV0.V1'][0])
    V[2] = np.abs(Solution['V3toV2.V2'][0])
    V[3] = np.abs(Solution['V3toV2.V3'][0]) 

    Vf = abs(Solution['V3toV2.V3'][2])
    Vr = abs(Solution['V3toV2.V3'][3])

    Gamma = np.abs(Vr/Vf)

    Vs = (V**2)/(Vf**2) - 1

    udenom = 2*(C[Bi]-C[Ci])*(CABBC - SABBC)
    vdenom = 2*(S[Bi]-S[Ci])*(SABBC - CABBC)

    u = ((Vs[Ai] - Vs[Bi]) - (Vs[Bi] - Vs[Ci])*SABBC)/udenom
    v = ((Vs[Ai] - Vs[Bi]) - (Vs[Bi] - Vs[Ci])*CABBC)/vdenom

    EpsB = 2*v/((1+u)**2 + v**2)
    EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))

    x = np.array([EpsG,EpsB])
    x = x/np.linalg.norm(x)
    EpsG = x[0]
    EpsB = x[1]

    return EpsG,EpsB


def ModAlgo4V(Solution,FREQ,ProbeIndexA=None,ProbeIndexB=None,ProbeIndexC=None):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([2.35,1.69,0.895,0.235])
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

    Vf = Solution['V3toV2.V3'][2]
    Vr = Solution['V3toV2.V3'][3]

    Gamma = Vr/Vf

    Vf = np.abs(Solution['V3toV2.V3'][2])

    Vs = (V**2)/(Vf**2)

    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

    EpsB = 2*v/((1+u)**2 + v**2)
    EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))

    beta = -1*np.angle(Gamma) #yolo
    x = np.array([np.cos(beta)*EpsG-np.sin(beta)*EpsB,np.sin(beta)*EpsG + np.cos(beta)*EpsB])
    x = x/np.linalg.norm(x)
    EpsG = x[0]
    EpsB = x[1]

    return EpsG,EpsB

def Algo4V(Solution,FREQ,ProbeIndexA=None,ProbeIndexB=None,ProbeIndexC=None):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    # Get voltages at positions in Measurepoints:
    MeasurePoints = np.array([2.35,1.69,0.895,0.235])
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


    Vs = (V**2)/(Vf**2)

    u = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigS)/((C[0] - C[1]) - (C[2] - C[3])*BigS)
    v = (1/2)*((Vs[0] - Vs[1]) - (Vs[2] - Vs[3])*BigC)/((S[0] - S[1]) - (S[2] - S[3])*BigC) 

    EpsB = 2*v/((1+u)**2 + v**2)
    EpsG = 1 - ((1-u**2-v**2)/((1+u)**2 + v**2))

    x = np.array([EpsG,EpsB])
    x = x/np.linalg.norm(x)
    EpsG = x[0]
    EpsB = x[1]

    return EpsG,EpsB

def Algo2V(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC=None):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    A = ProbeIndexA
    B = ProbeIndexB

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

    EpsB = (V[A] - Vf)*S[B] - (V[B] - Vf)*S[A]
    EpsG = (V[A] - Vf)*C[B] - (V[B] - Vf)*C[A]

    x = np.array([EpsG,EpsB])
    x = x/np.linalg.norm(x)
    EpsG = x[0]
    EpsB = x[1]

    return -1*EpsG,EpsB

def ModAlgo2V(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC=None):
    ######################################
    # Calculate new values of Cs and Cp  #
    ######################################

    A = ProbeIndexA
    B = ProbeIndexB

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

    Vf = Solution['V3toV2.V3'][2]
    Vr = Solution['V3toV2.V3'][3]

    Gamma = Vr/Vf

    Vf = np.abs(Vf)
    EpsB = (V[A] - Vf)*S[B] - (V[B] - Vf)*S[A]
    EpsG = (V[A] - Vf)*C[B] - (V[B] - Vf)*C[A]

    beta = -1*np.angle(Gamma) #yolo
    x = np.array([np.cos(beta)*EpsG-np.sin(beta)*EpsB,np.sin(beta)*EpsG + np.cos(beta)*EpsB])
    x = x/np.linalg.norm(x)
    EpsG = x[0]
    EpsB = x[1]

    return -1*EpsG,EpsB


