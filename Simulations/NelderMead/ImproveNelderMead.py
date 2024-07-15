import csv
import matplotlib.pyplot as plt
from matplotlib import cm,ticker
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import Circuit,TomasCircuit
from alive_progress import alive_bar
from scipy.stats import qmc

MHz_ = 1e6
global pF_
pF_ = 1e-12

global CsVals 
global CpVals 
global CaVals 

global NelderMeadCs 
global NelderMeadCp 

NelderMeadCs = []
NelderMeadCp = []

CsVals = np.linspace(800e-12,950e-12,20)#capacitor is defined in steps
CpVals = np.linspace(800e-12,950e-12,20)#capacitor is defined in steps
CaVals = np.linspace(235e-12,1000e-12,1)#capacitor is defined in steps
#leftmost = (np.abs(CsVals-25e-12)).argmin()
#rightmost = (np.abs(CsVals-400e-12)).argmin()
#CsVals = CsVals[leftmost:rightmost]

U = []
V = []


def MinimizeableFunction(CVals,CaVal,FREQ=25*MHz_):
    CpVal = CVals[0]
    CsVal = CVals[1]

    NelderMeadCs.append(CsVal)
    NelderMeadCp.append(CpVal)

    ct.set('Ca.Cs',CaVal)
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

    Vf = abs(Solution['V2toV3.V3'][2])
    Vr = abs(Solution['V2toV3.V3'][3])

    Gamma = Vr/Vf

    return Gamma

FREQ = 25*MHz_

FinalGammas = []

ct = TomasCircuit()

CaVal = 278*pF_

with alive_bar(len(CaVals)*len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    sampler = qmc.Halton(d=2, scramble=False)
    samples = sampler.random(n=1000)
    l_bounds = [7e-12,25e-12]
    u_bounds = [1000e-12,1000e-12]
    samples = qmc.scale(samples, l_bounds, u_bounds)
    for sample in samples:
        CpVal = sample[0]
        CsVal = sample[1]
        result = optimize.minimize(MinimizeableFunction, np.array([CpVal,CsVal]),args=(CaVal,FREQ),method='Nelder-Mead')
        CsVals = 

        FinalGammas.append(result.fun)

        bar()

    """
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):

            result = optimize.minimize(MinimizeableFunction, np.array([CpVal,CsVal]),args=(CaVal,FREQ),method='Nelder-Mead')

            FinalGammas.append(result.fun)

            bar()
    """
    X, Y = np.meshgrid(CpVals/pF_,CsVals/pF_) 
    """
    X has form [[CpVals],[CpVals],..]
    Y has form [[CsVals[0]],[CsVals[1]],..]
    I.e a MESHGRID is created from the x and y axes (as the name suggests)
    """
    Z = np.array(FinalGammas).reshape(len(CpVals),len(CpVals)) #put into shape with x length CpVals and y length CsVals
    levels = np.linspace(0.0, 1.0, 11)
    CS = plt.contourf(X, Y, Z,levels=levels, cmap=cm.coolwarm, extend='min')
    colorbar = plt.colorbar(CS)
    plt.ylabel("Cs (pF)")
    plt.title(f"Ca = {round(CaVal/pF_,2)}, Nelder mead matches i.f.o starting position")

plt.show()
