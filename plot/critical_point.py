import numpy as np
import matplotlib.pyplot as plt

b = 2
L = [32, 64, 128]
colors = ["b", "g", ""]
data = []

for l in L:

    filename = "data/critical_point_L_"+str(l)+"_K_-0.44.txt"
    data.append(np.loadtxt(filename, unpack=True, delimiter=","))
    
for i in range(len(L)):
    
    plt.plot(data[i][0], data[i][2], label="L = "+str(L[i])+", S = "+str(round(L[i]/b)))

plt.legend()
plt.show()
