import csv
import numpy as np
from pyRFtk import plotVSWs
from scipy.stats import qmc
import matplotlib.pyplot as plt
from matplotlib import cm,ticker
from scipy import spatial, optimize
from alive_progress import alive_bar
from TomasCircuit import Circuit,TomasCircuit
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler


#-------------------------------------------------------------------------------#
#                               Colors                                          #
#-------------------------------------------------------------------------------#

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


global MHz_
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


def MinimizeableFunction(CVals,FREQ=25*MHz_):
    CpVal = CVals[0]
    CsVal = CVals[1]
    CaVal = CVals[2]

    NelderMeadCs.append(CsVal)
    NelderMeadCp.append(CpVal)

    ct.set('Ca.Cs',CaVal)
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V0toV1','V2toV3','Ca','Cp','Cs'])

    Vf = abs(Solution['V2toV3.V3'][2])
    Vr = abs(Solution['V2toV3.V3'][3])

    Gamma = Vr/Vf

    NelderMeadGammas.append(Gamma)

    rICa = (1/np.sqrt(2))*np.abs(Solution['Ca.cap'][1])/138
    rICp = (1/np.sqrt(2))*np.abs(Solution['Cp.Cp'][1])/70
    rICs = (1/np.sqrt(2))*np.abs(Solution['Cs.cap'][1])/140
    IMaxR = max(np.array([rICa,rICp,rICs]))

    return Gamma + 0.01*IMaxR

FREQs = np.linspace(10,55,91)*MHz_

f = open('Data.csv','a')

writer = csv.writer(f)

sampler = qmc.Halton(d=3, scramble=False)
samples = sampler.random(n=100)
l_bounds = [7e-12,25e-12,35e-12]
u_bounds = [1000e-12,1000e-12,1000e-12]
samples = qmc.scale(samples, l_bounds, u_bounds)

with alive_bar(len(FREQs)*len(samples),title="Making colormap",length=20,bar="filling",spinner="dots_waves2") as bar:
    for FREQ in FREQs:

        Gammas = []

        ct = TomasCircuit()

        bounds = [(7e-12,1000e-12),(25e-12,1000e-12),(35e-12,1000e-12)]
        BestGamma = 1
        for sample in samples:
            CpVal = sample[0]
            CsVal = sample[1]
            CaVal = sample[2]
            result = optimize.minimize(MinimizeableFunction, np.array([CpVal,CsVal,CaVal]),args=(FREQ),bounds=bounds,method='Nelder-Mead')

            if BestGamma > result.fun:
                BestGamma = result.fun
                BestCpVal = result.x[0]
                BestCsVal = result.x[1]
                BestCaVal = result.x[2]
            bar()

        if BestGamma > 0.1:
            print(f"{bcolors.FAIL}No Gamma < 0.1 found...{bcolors.ENDC}")
        writer.writerow([FREQ/MHz_,BestCpVal,BestCsVal,BestCaVal])
        print([FREQ/MHz_,BestCpVal,BestCsVal,BestCaVal])
        print(f"Gamma: {BestGamma}")

f.close
