import csv
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy import spatial, optimize
from pyRFtk import rfCircuit, rfTRL, rfRLC, rfObject, rfCoupler
from pyRFtk import plotVSWs
from scipy.stats import qmc
from TomasCircuit import ArthurMod
from alive_progress import alive_bar
from Algos import *

"""
Finds Ca(f)
"""

plt.rcParams.update({'font.size': 20})

global MHz_
MHz_ = 1e6
global pF_
pF_ = 1e-12

###################
#Define capacitors#
###################
global mincs 
mincs = 25e-12
global maxcs 
maxcs = 1000e-12
global mincp 
mincp = 7e-12
global maxcp 
maxcp = 1000e-12
global minca
minca = 7e-12
global maxca
maxca = 1000e-12

global CaVals
CaVals = np.linspace(35e-12,1000e-12,8000)#capacitor is defined in steps
leftmost = (np.abs(CaVals-35e-12)).argmin()
rightmost = (np.abs(CaVals-200e-12)).argmin()#region of interest
CaVals = CaVals[leftmost:rightmost]
global CsVals
CsVals = np.linspace(mincs,maxcs,5000)
global CpVals
CpVals = np.linspace(mincp,maxcp,5000)
def InitialCa(FREQ):
    CaSugg = (2519*np.exp(-0.105*FREQ/MHz_) + 37 - 10)*pF_
    print(f"Starting at Ca={CaSugg/pF_}pF")
    return CaSugg
def AccMatch(ct,FREQ):
    l_bounds = [7e-12,25e-12]
    u_bounds = [1000e-12,1000e-12]
    options = {'maxiter':30}
    #Find Initial guess through matching algo
    CpVal,CsVal = Match(DirCoupler,ct,FREQ)
    #Find best of the best through Nelder mead
    result = optimize.minimize(MinimizeableFunction, np.array([CpVal,CsVal]),args=(FREQ),method='Nelder-Mead',options=options)
    Gamma = result.fun
    CVals = result.x
    CpVal,CsVal = CVals[0],CVals[1]
    return CpVal,CsVal
def CurrentEstimate(ct,CpVal,CsVal):
    #Define small square around point
    SquareLength = 4e-12

    leftmost = (np.abs(CpVals-(CpVal-SquareLength/2))).argmin()
    rightmost = (np.abs(CpVals-(CpVal+SquareLength/2))).argmin()
    CpValRegion = CpVals[leftmost:rightmost]

    leftmost = (np.abs(CsVals-(CsVal-SquareLength/2))).argmin()
    rightmost = (np.abs(CsVals-(CsVal+SquareLength/2))).argmin()
    CsValRegion = CsVals[leftmost:rightmost]
    
    points = 0
    Iratios = []

    for CsVal in CsValRegion:
        for CpVal in CpValRegion:
            CVals = [CpVal,CsVal]
            Gamma = MinimizeableFunction(CVals,FREQ)
            if Gamma < 0.1:
                Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*1.0*1E3)}, nodes=['V1toV0','V3toV2','Cs','Cp','Ca'])
                rICa = (1/np.sqrt(2))*np.abs(Solution['Ca.cap'][1])/(138*(1-Gamma))
                rICp = (1/np.sqrt(2))*np.abs(Solution['Cp.Cp'][1])/(70*(1-Gamma))
                rICs = (1/np.sqrt(2))*np.abs(Solution['Cs.cap'][1])/(140*(1-Gamma))
                Iratios.append(max(np.array([rICa,rICp,rICs])))
                points += 1
    try:
        MaxI = max(np.array(Iratios))
    except:
        MaxI = 400
    return points,MaxI

def CurrentWeightedArea(CaVal,FREQ,ct):
    CaVal = CaVal[0]
    ct.set('Ca.Cs',CaVal)
    ct.set('Cp.Cs',100*pF_)
    ct.set('Cs.Cs',100*pF_)

    #########################################
    #           Match the system            #
    #########################################
    CpVal,CsVal = AccMatch(ct,FREQ)
    #########################################
    #Get current and estimate of domain size#
    #########################################
    points,MaxI = CurrentEstimate(ct,CpVal,CsVal)
    if points > 6: 
        return MaxI
    else:
        return 400


def MinimizeableFunction(CVals,FREQ=25*MHz_):
    CpVal = CVals[0]
    CsVal = CVals[1]
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V1toV0','V3toV2'])

    Vf = abs(Solution['V3toV2.V3'][2])
    Vr = abs(Solution['V3toV2.V3'][3])

    Gamma = Vr/Vf

    return Gamma

def Match(algorithm,ct,FREQ):
    matched = False
    j=0
    while not matched and j < 4:
        CsVal = np.random.choice(CsVals)
        CpVal = np.random.choice(CpVals)
        if j==1:
            CpVal = 200*pF_
            CsVal = 40*pF_
        i=0
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

            if abs(Vr)/abs(Vf) <= 0.05:
                matched = True
                return CpVal,CsVal
            ProbeIndexA,ProbeIndexB,ProbeIndexC = 0,1,3
            EpsG,EpsB = algorithm(Solution,FREQ,ProbeIndexA,ProbeIndexB,ProbeIndexC)
            
            SPFactor = 5*pF_
            if abs(Gamma) < 0.1: SPFactor = 1*pF_
            CsVal += SPFactor*EpsG
            CpVal += SPFactor*EpsB

            i += 1
        j+=1
    return 100e-12,100e-12




FREQs = [50*MHz_]

# CSV file
f = open(f'Ca_Hf.csv', 'w')
writer = csv.writer(f)
writer.writerow(["Frequency (MHz)","Ideal CaVal (pf)","Middle CsVal (pf)","Middle CpVal (pf)","MaxI"])

with alive_bar(len(FREQs),title="Finding #Matched/#Total",length=20,bar="filling",spinner="dots_waves2") as bar:
    for FREQ in FREQs:
        ct = ArthurMod()
        #options = {'maxiter':10,'xatol':5*pF_}
        options = {}
        cabounds = optimize.Bounds(minca,maxca)
        result = optimize.minimize(CurrentWeightedArea,InitialCa(FREQ),args=(FREQ,ct),bounds=cabounds,method='Nelder-Mead',options=options)
        CaVal = result.x[0]
        CpVal,CsVal = AccMatch(ct,FREQ)
        if result.fun < 10:
            _,MaxI = CurrentEstimate(ct,CpVal,CsVal)
            writer.writerow([FREQ/MHz_,CaVal/pF_,CsVal/pF_,CpVal/pF_,MaxI])
            print(f"f={FREQ},Ca:{CaVal},Cs:{CsVal},Cp:{CpVal}")
        else:
            print(f"f={FREQ} Failed, final Gamma: {MinimizeableFunction([CpVal,CsVal],FREQ)}")
        bar()
f.close()
