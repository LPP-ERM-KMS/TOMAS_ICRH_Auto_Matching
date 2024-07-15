import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from TomasCircuit import TomasCircuit
from alive_progress import alive_bar

MHz_ = 1e6
pF_ = 1e-12

def findClosestValue(givenList, target):
    def difference(givenList):
        return abs(givenList - target)
    result = min(givenList, key=difference)
    return result

def ComputeGamma(CVals,FREQ=25*MHz_):
    CaVal = CVals[0]
    CpVal = CVals[1]
    CsVal = CVals[2]
    CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
    CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps
    CaVals = np.linspace(35e-12,1000e-12,8000)#capacitor is defined in steps

    CpVal = CpVals[(np.abs(CpVals-CpVal)).argmin()]
    CsVal = CsVals[(np.abs(CsVals-CsVal)).argmin()]
    CaVal = CaVals[(np.abs(CaVals-CaVal)).argmin()]

    ct.set('Ca.Cs',CaVal)
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3'])

    Vf = abs(Solution['V2toV3.V3'][2])
    Vr = abs(Solution['V2toV3.V3'][3])

    Gamma = np.abs(Vr/Vf)
    return Gamma



FREQ = 25.0*MHz_
ct = TomasCircuit()

print(optimize.minimize(ComputeGamma, np.array([200*pF_,200*pF_,200*pF_]),method='Nelder-Mead'))
