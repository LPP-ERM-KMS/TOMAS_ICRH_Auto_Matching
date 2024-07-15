import csv
import matplotlib.pyplot as plt
from matplotlib import cm,ticker
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import Circuit,TomasCircuit
from alive_progress import alive_bar

MHz_ = 1e6
global pF_
pF_ = 1e-12

global NelderMeadCs 
global NelderMeadCp 

NelderMeadCs = []
NelderMeadCp = []
NelderMeadGammas = []

CsVals = np.linspace(25e-12,1000e-12,500)#capacitor is defined in steps
leftmost = (np.abs(CsVals-800e-12)).argmin()
rightmost = (np.abs(CsVals-950e-12)).argmin()
CsVals = CsVals[leftmost:rightmost]

CpVals = np.linspace(7e-12,1000e-12,500)#capacitor is defined in steps
leftmost = (np.abs(CpVals-800e-12)).argmin()
rightmost = (np.abs(CpVals-950e-12)).argmin()
CpVals = CpVals[leftmost:rightmost]

CaVals = np.linspace(235e-12,1000e-12,1)#capacitor is defined in steps

U = []
V = []


def ComputeGamma(CVals,FREQ=25*MHz_):
    CpVal = CVals[0]
    CsVal = CVals[1]
    CaVal = CVals[2]

    #CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
    #CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps

    #CpVal = CpVals[(np.abs(CpVals-CpVal)).argmin()]
    #CsVal = CsVals[(np.abs(CsVals-CsVal)).argmin()]
        
    NelderMeadCs.append(CsVal)
    NelderMeadCp.append(CpVal)

    ct.set('Ca.Cs',CaVal)
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

    Vf = abs(Solution['V2toV3.V3'][2])
    Vr = abs(Solution['V2toV3.V3'][3])

    Gamma = Vr/Vf

    NelderMeadGammas.append(Gamma)
    return Gamma

FREQ = 25*MHz_

Gammas = []

ct = TomasCircuit()
#ct = Circuit(WhichStrap=3)

CaVal = 276.5*pF_
"""
with alive_bar(len(CaVals)*len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="halloween",spinner="dots_waves2") as bar:
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):
            ct.set('Ca.Cs',CaVal)
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)

            
            Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])
            Vf = abs(Solution['V2toV3.V3'][2])
            Vr = abs(Solution['V2toV3.V3'][3])

            Gammas.append(Vr/Vf)

            bar()

    X, Y = np.meshgrid(CpVals/pF_,CsVals/pF_) 
"""
#    X has form [[CpVals],[CpVals],..]
#    Y has form [[CsVals[0]],[CsVals[1]],..]
#    I.e a MESHGRID is created from the x and y axes (as the name suggests)
"""
    Z = np.array(Gammas).reshape(len(CsVals),len(CpVals)) #put into shape with x length CpVals and y length CsVals
    levels = np.linspace(0.0, 1.0, 11)
    CS = plt.contourf(X, Y, Z,levels=levels, cmap=cm.coolwarm, extend='min')
    colorbar = plt.colorbar(CS)
    plt.ylabel("Cs (pF)")
    plt.xlabel("Cp (pF)")
    plt.title(f"Ca = {round(CaVal/pF_,2)}")
"""
bounds = [(25e-12,1000e-12),(7e-12,1000e-12),(35e-12,1000e-12)]
result = optimize.minimize(ComputeGamma, np.array([200*pF_,200*pF_,200*pF_]),args=(FREQ),bounds=bounds,method='Nelder-Mead',options={'return_all':True})
##plt.scatter(np.array(result.allvecs)[:,0]/pF_,np.array(result.allvecs)[:,1]/pF_)
##plt.plot(np.array(result.allvecs)[:,0]/pF_,np.array(result.allvecs)[:,1]/pF_)
##plt.show()
print(result.fun)
print(result.x)
"""
plt.show()
Gammas = []
for vector in result.allvecs:
    CVals = np.zeros(2)
    CVals[0] = vector[0]
    CVals[1] = vector[1]
    Gammas.append(ComputeGamma(CVals,CaVal,FREQ=FREQ))

plt.plot(Gammas)
plt.show()
"""

