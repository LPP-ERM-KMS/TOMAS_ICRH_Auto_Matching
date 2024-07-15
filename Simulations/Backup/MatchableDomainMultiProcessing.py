import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
from alive_progress import alive_bar
import concurrent.futures

def wave(Gamma,beta,x):
    return np.abs(np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta))

def GammaErrors(Gamma,f,Probe1,Probe2,Vf=1):
    c = 3e8
    VoltMeterPositions = np.array([0.235,0.895,1.69,2.35,0])
    x = np.linspace(0,2.59,1000)
    Y = (1-Gamma)/(1+Gamma)
    lam = f/c
    beta = 2*np.pi*lam
    V = np.zeros(len(VoltMeterPositions))
    for k,VoltMeterPosition in enumerate(VoltMeterPositions):
        V[k] = wave(Gamma,beta,x)[np.abs(x-VoltMeterPosition).argmin()]
    S = np.sin(2*beta*VoltMeterPositions) #array
    C = np.cos(2*beta*VoltMeterPositions) #array
    IYCorr = (V[Probe1] - Vf)*S[Probe2] - (V[Probe2] - Vf)*S[Probe1]
    RYCorr = (V[Probe1] - Vf)*C[Probe2] - (V[Probe2] - Vf)*C[Probe1]
    Y_n = Y + RYCorr + 1j*IYCorr
    NewGamma = (1 - Y_n)/(1 + Y_n)
    GammaRealError = NewGamma.real - Gamma.real
    GammaImagError = NewGamma.imag - Gamma.imag
    return GammaRealError,GammaImagError

def Matchable(InitialPoint,f,Probe1,Probe2,ScaleFactor=0.025):
    Gamma = InitialPoint

    step = 0
    Matched = False
    while not Matched:        
        RCorr,CCorr = GammaErrors(Gamma,f,Probe1,Probe2)
        Gamma += ScaleFactor*(RCorr + 1j*CCorr)
        step += 1
        if step > 2000 or np.abs(Gamma) > 1:
            break
        if np.abs(Gamma) < 0.2:
            Matched = True

    return Matched


fs = np.linspace(25e6,55e6,30)
MatchRatio = np.zeros((len(fs),4,4,1))
Resolution = 100

with alive_bar(len(fs)*Resolution*Resolution,title='Calculating',length=20,bar='filling',spinner='dots_waves2') as bar:
    for fiter,f in enumerate(fs):
        for i in np.linspace(-1,1,Resolution):
            for j in np.linspace(0,1,int(Resolution/2)):
                bar()
                if i**2 + j**2 > 1:
                    continue

                #Multithreading
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    t01p = executor.submit(Matchable,i+1j*j,f,0,1)
                    t02p = executor.submit(Matchable,i+1j*j,f,0,2)
                    t03p = executor.submit(Matchable,i+1j*j,f,0,3)
                    t12p = executor.submit(Matchable,i+1j*j,f,1,2)
                    t13p = executor.submit(Matchable,i+1j*j,f,1,3)
                    t23p = executor.submit(Matchable,i+1j*j,f,2,3)

                    t01n = executor.submit(Matchable,i-1j*j,f,0,1)
                    t02n = executor.submit(Matchable,i-1j*j,f,0,2)
                    t03n = executor.submit(Matchable,i-1j*j,f,0,3)
                    t12n = executor.submit(Matchable,i-1j*j,f,1,2)
                    t13n = executor.submit(Matchable,i-1j*j,f,1,3)
                    t23n = executor.submit(Matchable,i-1j*j,f,2,3)

                    if t01p.result():
                        MatchRatio[fiter][0][1] = MatchRatio[fiter][0][1] + 1/(Resolution*Resolution)
                    if t02p.result():
                        MatchRatio[fiter][0][2] = MatchRatio[fiter][0][2] + 1/(Resolution*Resolution)
                    if t03p.result():
                        MatchRatio[fiter][0][3] = MatchRatio[fiter][0][3] + 1/(Resolution*Resolution)
                    if t12p.result():
                        MatchRatio[fiter][1][2] = MatchRatio[fiter][1][2] + 1/(Resolution*Resolution)
                    if t13p.result():
                        MatchRatio[fiter][1][3] = MatchRatio[fiter][1][3] + 1/(Resolution*Resolution)
                    if t23p.result():
                        MatchRatio[fiter][2][3] = MatchRatio[fiter][2][3] + 1/(Resolution*Resolution)

                    if t01n.result():
                        MatchRatio[fiter][0][1] = MatchRatio[fiter][0][1] + 1/(Resolution*Resolution)
                    if t02n.result():
                        MatchRatio[fiter][0][2] = MatchRatio[fiter][0][2] + 1/(Resolution*Resolution)
                    if t03n.result():
                        MatchRatio[fiter][0][3] = MatchRatio[fiter][0][3] + 1/(Resolution*Resolution)
                    if t12n.result():
                        MatchRatio[fiter][1][2] = MatchRatio[fiter][1][2] + 1/(Resolution*Resolution)
                    if t13n.result():
                        MatchRatio[fiter][1][3] = MatchRatio[fiter][1][3] + 1/(Resolution*Resolution)
                    if t23n.result():
                        MatchRatio[fiter][2][3] = MatchRatio[fiter][2][3] + 1/(Resolution*Resolution)

                bar()

for k in range(0,4):
    for l in range(0,4):
        if k<l:
            plt.plot(fs/1e6,MatchRatio[:,k,l],label=f"Probe {k} and Probe {l}")

plt.xlabel("Frequency (MHz)")
plt.ylabel("Matchable surface ratio")
plt.legend()
plt.show()
