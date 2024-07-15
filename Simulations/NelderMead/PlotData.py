import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import genfromtxt

"""
Explores the domain and looks at the size of the blue blob and the maximal current running through the capacitors inside the blob
"""

plt.rcParams.update({'font.size': 20})

MHz_=1e6
pF_=1e-12
FREQ = 25*MHz_

MyData = genfromtxt('data.csv', delimiter=',')
#row = "CaVal","points","MaxI" 
print(MyData)
CaVals = MyData[:,0]
TotalPoints = MyData[:,1]
MaxI = MyData[:,2]

COLOR_POINTS = "#69b3a2"
COLOR_CURRENT = "#3399e6"

fig, ax1 = plt.subplots(figsize=(8, 8))
ax2 = ax1.twinx()

ax1.plot(CaVals/pF_, TotalPoints, color=COLOR_POINTS, lw=3)
ax2.plot(CaVals/pF_, MaxI, color=COLOR_CURRENT, lw=4)

ax1.set_xlabel("Ca (pF)")
ax1.set_ylabel("Cs Cp combinations leading to $\Gamma<0.1$", color=COLOR_POINTS)
ax1.tick_params(axis="y", labelcolor=COLOR_POINTS)

ax2.set_ylabel("Ratio of maximum observed RMS current to the respective capacitor limit", color=COLOR_CURRENT)
ax2.tick_params(axis="y", labelcolor=COLOR_CURRENT)

plt.title(f"TOMAS ICRH circuit at 1kW {FREQ/MHz_}MHz")
plt.savefig(f"{FREQ/MHz_}MHzCaScan.pdf")
plt.show()
