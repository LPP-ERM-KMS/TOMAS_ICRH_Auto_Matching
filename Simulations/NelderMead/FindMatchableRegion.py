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

"""
Explores the domain and looks at the size of the blue blob and the maximal current running through the capacitors inside the blob
"""

plt.rcParams.update({'font.size': 20})

MHz_ = 1e6
pF_ = 1e-12

def MinimizeableFunction(CVals,CaVal,FREQ=25*MHz_):
    CpVal = CVals[0]
    CsVal = CVals[1]

    ct.set('Ca.Cs',CaVal)
    ct.set('Cs.Cs',CsVal)
    ct.set('Cp.Cs',CpVal)
    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*5*1E3)}, nodes=['V1toV0','V3toV2'])

    Vf = abs(Solution['V3toV2.V3'][2])
    Vr = abs(Solution['V3toV2.V3'][3])

    Gamma = Vr/Vf

    return Gamma


CaVals = np.linspace(35e-12,1000e-12,8)#capacitor is defined in steps
leftmost = (np.abs(CaVals-35e-12)).argmin()
rightmost = (np.abs(CaVals-800e-12)).argmin()
CaVals = CaVals[leftmost:rightmost]

FREQ = 25*MHz_
NPoints = 120 #culprit, if not high enough, won't find the ideal region random points

TotalPoints = []
MaxI = []

f = open('data50MHz.csv', 'w')
# create the csv writer
writer = csv.writer(f)

writer.writerow(["CaVal","points","MaxI"])

ct = ArthurMod()

with alive_bar(len(CaVals)*NPoints + len(CaVals)*409*401,title="Finding #Matched/#Total",length=20,bar="filling",spinner="dots_waves2") as bar: #410,402
    for CaVal in CaVals:

        CsVals = np.linspace(25e-12,1000e-12,5000)#capacitor is defined in steps
        CpVals = np.linspace(7e-12,1000e-12,5000)#capacitor is defined in steps

        l_bounds = [7e-12,25e-12]
        u_bounds = [1000e-12,1000e-12]

        sampler = qmc.Halton(d=2, scramble=False)
        samples = sampler.random(n=NPoints)
        samples = qmc.scale(samples, l_bounds, u_bounds)

        options = {'maxiter':30}

        FinalGammas = []
        FinalCaps = []

        for i,sample in enumerate(samples):
            CpVal = sample[0]
            CsVal = sample[1]
            result = optimize.minimize(MinimizeableFunction, np.array([CpVal,CsVal]),args=(CaVal,FREQ),method='Nelder-Mead',options=options)
            FinalGammas.append(result.fun)
            FinalCaps.append(result.x)
            bar()

        BestCaps = FinalCaps[np.where(np.array(FinalGammas)==min(np.array(FinalGammas)))[0][0]]

        print(min(np.array(FinalGammas)))
        print(BestCaps)
        print(CaVal)

        SquareLength = 80e-12 #ideally 80e-12

        leftmost = (np.abs(CpVals-(BestCaps[0]-SquareLength/2))).argmin()
        rightmost = (np.abs(CpVals-(BestCaps[0]+SquareLength/2))).argmin()
        CpValRegion = CpVals[leftmost:rightmost]

        leftmost = (np.abs(CsVals-(BestCaps[1]-SquareLength/2))).argmin()
        rightmost = (np.abs(CsVals-(BestCaps[1]+SquareLength/2))).argmin()
        CsValRegion = CsVals[leftmost:rightmost]
        
        points = 0
        Iratios = []

        for CsVal in CsValRegion:
            for CpVal in CpValRegion:
                CVals = [CpVal,CsVal]
                Gamma = MinimizeableFunction(CVals,CaVal,FREQ)
                if Gamma < 0.1:
                    Solution = ct.Solution(f=FREQ, E={'Source': np.sqrt(2*50*1.0*1E3)}, nodes=['V1toV0','V3toV2','Cs','Cp','Ca'])
                    rICa = (1/np.sqrt(2))*np.abs(Solution['Ca.cap'][1])/(138*(1-Gamma))
                    rICp = (1/np.sqrt(2))*np.abs(Solution['Cp.Cp'][1])/(70*(1-Gamma))
                    rICs = (1/np.sqrt(2))*np.abs(Solution['Cs.cap'][1])/(140*(1-Gamma))
                    Iratios.append(max(np.array([rICa,rICp,rICs])))
                    points += 1
                bar()
        TotalPoints.append(points)
        try:
            MaxI.append(max(np.array(Iratios)))
        except:
            print(f"no matchable region found for Ca={CaVal/pF_}pF")
            MaxI.append(0)
        writer.writerow([CaVal,points,MaxI[-1]])
        
f.close()
COLOR_POINTS = "#69b3a2"
COLOR_CURRENT = "#3399e6"

fig, ax1 = plt.subplots(figsize=(8, 8))
ax2 = ax1.twinx()

ax1.plot(CaVals/pF_, TotalPoints, color=COLOR_POINTS, lw=3)
ax2.plot(CaVals/pF_, MaxI, color=COLOR_CURRENT, lw=4)

ax1.set_xlabel("Ca (pF)")
ax1.set_ylabel("Cs Cp combinations leading to $\Gamma<0.1$", color=COLOR_POINTS)
ax1.tick_params(axis="y", labelcolor=COLOR_POINTS)

ax2.set_ylabel("Maximum RMS current ratio to maximum of a capacitor in the matched region", color=COLOR_CURRENT)
ax2.tick_params(axis="y", labelcolor=COLOR_CURRENT)

plt.title(f"TOMAS ICRH circuit at 1kW {FREQ/MHz_}MHz")
plt.savefig(f"{FREQ/MHz_}MHzCaScan.pdf")
plt.show()
