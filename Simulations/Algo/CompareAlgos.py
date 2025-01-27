import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import ArthurMod as TomasCircuit
from alive_progress import alive_bar
from Algos import *
from mpl_smithchart import SmithAxes


MHz_ = 1e6
pF_ = 1e-12

resolution = 12
maxcs = 1000e-12
mincs = 7e-12
mincp = 35e-12
maxcp = 1000e-12

FREQ = 25*MHz_

CsVals = np.linspace(mincs,maxcs,resolution)#capacitor is defined in steps
#leftmost = (np.abs(CsVals-120e-12)).argmin()
#rightmost = (np.abs(CsVals-170e-12)).argmin()
#CsVals = CsVals[leftmost:rightmost]

CpVals = np.linspace(mincp,maxcp,resolution)#capacitor is defined in steps
#leftmost = (np.abs(CpVals-100e-12)).argmin()
#rightmost = (np.abs(CpVals-750e-12)).argmin()
#CpVals = CpVals[leftmost:rightmost]

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
def MatchAble(algorithm,ct,CpVal,CsVal,FREQ,ProbeIndexA=None,ProbeIndexB=None,ProbeIndexC=None):
    i = 0
    matched = False
    while ((i < 10000) and not matched):
        
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

        if abs(Vr)/abs(Vf) <= 0.1:
            matched = True
            print("Match")

        EpsG,EpsB = algorithm(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC)
        
        CsVal += EpsG*5*pF_
        CpVal += EpsB*5*pF_

        i += 1

    return matched


Gammas = []

ct = TomasCircuit()

CaVal = 150*pF_
ct.set('Ca.Cs',CaVal)

num = 0
tot = 0

with alive_bar(len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)
            matched = MatchAble(DirCoupler,ct,CpVal,CsVal,FREQ,0,2,3)
            tot+=1
            if matched:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='green')
                num+=1
            else:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='red')
            bar()

plt.ylabel("Cs (pF)")
plt.xlabel("Cp (pF)")
plt.title(f"Directional Coupler Algorithm at {FREQ/MHz_}MHz, matchable Fraction = {num/tot}")
plt.savefig(f"DirCoupler_{FREQ/MHz_}.pdf")
plt.close()

"""
num = 0
tot = 0

with alive_bar(len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)
            matched = MatchAble(Algo3V,ct,CpVal,CsVal,FREQ,0,2,3)
            tot+=1
            if matched:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='green')
                num+=1
            else:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='red')
            bar()

plt.ylabel("Cs (pF)")
plt.xlabel("Cp (pF)")
plt.title(f"3V Algorithm at {FREQ/MHz_}MHz. Matchable Fraction = {num/tot}")
plt.savefig(f"3VMatchability_{FREQ/MHz_}.pdf")
plt.close()

num = 0
tot = 0

with alive_bar(len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)
            matched = MatchAble(Algo4V,ct,CpVal,CsVal,FREQ,0,2,3)
            tot+=1
            if matched:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='green')
                num+=1
            else:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='red')
            bar()

plt.ylabel("Cs (pF)")
plt.xlabel("Cp (pF)")
plt.title(f"4V Algorithm at {FREQ/MHz_}MHz. Matchable Fraction = {num/tot}")
plt.savefig(f"4VMatchability_{FREQ/MHz_}.pdf")
plt.close()


num = 0
tot = 0

with alive_bar(len(CsVals)*len(CpVals),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for j,CsVal in enumerate(CsVals):
        for k,CpVal in enumerate(CpVals):
            ct.set('Cs.Cs',CsVal)
            ct.set('Cp.Cs',CpVal)
            matched = MatchAble(Algo2V,ct,CpVal,CsVal,FREQ,0,3,None)
            tot+=1
            if matched:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='green')
                num+=1
            else:
                plt.scatter(CpVal/pF_,CsVal/pF_,color='red')
            bar()

plt.ylabel("Cs (pF)")
plt.xlabel("Cp (pF)")
plt.title(f"2V Algorithm at {FREQ/MHz_}MHz. Matchable Fraction = {num/tot}")
plt.savefig(f"2VMatchability_{FREQ/MHz_}.pdf")
"""
