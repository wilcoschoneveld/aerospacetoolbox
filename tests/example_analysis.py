import numpy as np
import matplotlib.pyplot as plt
from aerospacetoolbox import *

#import the windtunnel distribution data (from the HSWTT reader, Table 2)
dist = np.matrix("""
[19.8 30.112 1.843; 24.8 27.325 1.673;
29.8 24.832 1.520; 34.8 22.642 1.386;
39.8 20.763 1.271; 44.8 19.204 1.175;
56.8 16.819 1.029; 59.8 16.532 1.012;
62.8 16.372 1.002; 74.8 16.958 1.038;
79.8 17.664 1.081; 84.8 18.553 1.136;
96.8 21.067 1.290; 99.8 21.698 1.328;
102.8 22.306 1.365; 114.8 24.488 1.499;
119.8 25.275 1.547; 124.8 25.992 1.591;
129.8 26.638 1.631; 134.8 27.219 1.666;
139.8 27.734 1.698; 144.8 28.188 1.725;
149.8 28.584 1.750; 154.8 28.924 1.770;
159.8 29.212 1.788; 164.8 29.450 1.803;
176.8 29.832 1.826; 179.8 29.890 1.830;
182.8 29.935 1.832; 194.8 30.000 1.836]""")

#import the windtunnel experimental data (obtained during the HSWTT)
data1 = np.loadtxt("hswtdata_example.txt", skiprows=5)

#use the distribution data to calculate the theoretical isentropic flow
[m1, t1, p1, r1, a1] = flowisentropic(sub=dist[:,2])
[m2, t2, p2, r2, a2] = flowisentropic(sup=dist[:,2])

#Stich the subsonic and supersonic solutions together
psubsup = np.vstack((p1[:9],p2[9:]))

#plot the theoretical and experimental data
plt.plot(dist[:,0], psubsup, 'k--', label="Theory")
plt.plot(data1[:,0], data1[:,1], 'k^-', label="Experiment")
plt.axis([19.8, 194.8, 0, 1])
plt.legend(loc="upper right")
plt.xlabel("location $x$ [mm]")
plt.ylabel(r"pressure ratio $\frac{p}{p_0}$ [-]")
plt.show()

