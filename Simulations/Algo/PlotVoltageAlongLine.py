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

25.0,5.385007663975487e-10,1.39070071409994e-10,2.164126717225259e-10

FREQ = 25*MHz_

Gammas = []

ct = TomasCircuit()

CaVal = 2.164126717225259e-10
CsVal = 1.39070071409994e-10
CpVal = 5.385007663975487e-10

ct.set('Ca.Cs',CaVal)
ct.set('Cs.Cs',CsVal)
ct.set('Cp.Cs',CpVal)

Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

maxV, where, VSWs = ct.maxV(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)})
plotVSWs(VSWs)

plt.show()
