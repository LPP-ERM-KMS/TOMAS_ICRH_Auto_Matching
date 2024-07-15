import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod 
from alive_progress import alive_bar
from Algos import *


MHz_ = 1e6
pF_ = 1e-12

resolution = 20
maxcp = 500e-12
maxcs = 200e-12
mincp = 50e-12
mincs = 25e-12

CsVals = np.linspace(mincs,maxcs,resolution)#capacitor is defined in steps
leftmost = (np.abs(CsVals-0e-12)).argmin()
rightmost = (np.abs(CsVals-750e-12)).argmin()
CsVals = CsVals[leftmost:rightmost]

CpVals = np.linspace(mincp,maxcp,resolution)#capacitor is defined in steps
leftmost = (np.abs(CpVals-0e-12)).argmin()
rightmost = (np.abs(CpVals-750e-12)).argmin()
CpVals = CpVals[leftmost:rightmost]

U = []
V = []
def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )
def MatchPath(algorithm,ct,CpVal,CsVal,FREQ,ProbeIndexA=None,ProbeIndexB=None,ProbeIndexC=None,phasefactor=None):
    global admittances
    global gammas
    admittances = []
    gammas = []
    i = 0
    matched = False
    Path = []
    while ((i < 2000) and not matched):
        
        Path.append(np.array([CpVal,CsVal])/pF_)

        if CsVal > maxcs or CsVal < mincs:
            break
        if CpVal > maxcp or CpVal < mincp:
            break

        ct.set('Cs.Cs',CsVal)
        ct.set('Cp.Cs',CpVal)

        Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V1toV0','V3toV2'])

        Vf = Solution['V3toV2.V3'][2]
        Vr = Solution['V3toV2.V3'][3]
        Gamma = Vr/Vf
        gammas.append(Gamma)

        admittances.append((1-Gamma)/(1+Gamma))

        if abs(Vr/Vf) <= 0.05:
            matched = True

        if phasefactor:
            EpsG,EpsB = algorithm(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC,phasefactor)
        else:
            EpsG,EpsB = algorithm(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC)

        CsVal += EpsG*1*pF_
        CpVal += EpsB*1*pF_

        i += 1

    Path.append(np.array([CpVal,CsVal])/pF_)
    
    return np.array(Path)

FREQ = 30*MHz_

Gammas = []

ct = ArthurMod()

CaVal = 150*pF_
ct.set('Ca.Cs',CaVal)
Matched = []

with alive_bar(len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)
            Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V3toV2'])
            Vf = Solution['V3toV2.V3'][2]
            Vr = Solution['V3toV2.V3'][3]
            Gammas.append(Vr/Vf)
            bar()

X, Y = np.meshgrid(CpVals/pF_,CsVals/pF_) 
Gammas = np.array(Gammas)
Z = np.abs(Gammas).reshape(len(CsVals),len(CpVals)) #put into shape with x length CsVals and y length CpVals
levels = np.linspace(0.0, 1.0, 11)
CS = plt.contourf(X, Y, Z, levels=levels, cmap=cm.coolwarm, extend='min')
colorbar = plt.colorbar(CS)
plt.ylabel("Cs (pF)")
plt.xlabel("Cp (pF)")
plt.title(f"Ca = {round(CaVal/pF_,2)}, Freq= {FREQ/MHz_}MHz")

CsInit = 150*pF_
CpInit = 100*pF_

Path = MatchPath(Algo2V,ct,CsInit,CpInit,FREQ,0,3)
PathPlot = plt.plot(Path[:, 0], Path[:, 1],label="TEXTOR",color="yellow")[0]
print(f"2V Steps: {len(Path)}")
add_arrow(PathPlot)

Path = MatchPath(ModAlgo3V,ct,CsInit,CpInit,FREQ,0,2,3,-1)
PathPlot = plt.plot(Path[:, 0], Path[:, 1],label="3V Algorithm",color="green")[0]
print(f"Mod3V Steps: {len(Path)}")
add_arrow(PathPlot)

Path = MatchPath(LazyModAlgo,ct,CsInit,CpInit,FREQ)
PathPlot = plt.plot(Path[:, 0], Path[:, 1],label="Dir coupler Algorithm",color="purple")[0]
print(f"Dir coupler Steps: {len(Path)}")
add_arrow(PathPlot)

Path = MatchPath(ModAlgo4V,ct,CsInit,CpInit,FREQ,0,2,3)
PathPlot = plt.plot(Path[:, 0], Path[:, 1],label="4V Algorithm",color="blue")[0]
print(f"Mod4V Steps: {len(Path)}")
add_arrow(PathPlot)

plt.legend()
plt.show()
