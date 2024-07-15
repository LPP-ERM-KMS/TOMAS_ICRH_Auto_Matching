import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
from alive_progress import alive_bar
from numba import jit, cuda

@jit(target_backend='cuda')
def wave(Gamma,beta,x):
    return np.abs(np.exp(1j*x*beta) + Gamma*np.exp(-1j*x*beta))

@jit(target_backend='cuda')
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

@jit(target_backend='cuda')
def Matchable(InitialPoint,f,Probe1,Probe2,ScaleFactor=0.050):
    Gamma = InitialPoint

    step = 0
    Matched = False
    while not Matched:        
        RCorr,CCorr = GammaErrors(Gamma,f,Probe1,Probe2)
        Gamma += ScaleFactor*(RCorr + 1j*CCorr)
        step += 1
        if step > 1000 or np.abs(Gamma) > 1:
            break
        if np.abs(Gamma) < 0.2:
            Matched = True

    return Matched


fs = np.linspace(25e6,55e6,2)
MatchRatio = np.zeros((len(fs),4,4,1))
Resolution = 50

with alive_bar(len(fs)*Resolution*Resolution,title='Calculating',length=20,bar='filling',spinner='dots_waves2') as bar:
    for fiter,f in enumerate(fs):
        for i in np.linspace(-1,1,Resolution):
            for j in np.linspace(-1,1,Resolution):
                bar()
                if i**2 + j**2 > 1:
                    continue
                #Now loop over possible probe combinations:
                for k in range(0,4):
                    for l in range(0,4):
                        if k<l:
                            Matchable_ = Matchable(i + 1j*j,f,k,l)
                            if Matchable_:
                                MatchRatio[fiter][k][l] = MatchRatio[fiter][k][l] + 1/(Resolution*Resolution)
                            else:
                                MatchRatio[fiter][k][l] = MatchRatio[fiter][k][l] - 1/(Resolution*Resolution)

for k in range(0,4):
    for l in range(0,4):
        if k<l:
            plt.plot(fs,MatchRatio[:,k,l],label=f"Probe {i} and Probe {j}")

plt.xlabel("Frequency (MHz)")
plt.ylabel("Matchable surface ratio")
plt.legend()
plt.show()
