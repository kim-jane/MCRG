import numpy as np
import matplotlib.pyplot as plt

b = 2
L = [16, 32, 64, 128]
colors = ["indianred", "orange", "forestgreen", "royalblue"]
data = []

plt.figure(figsize=(8,6))

Kc = -np.log(1+np.sqrt(2))/2
plt.hlines(Kc, xmin=-10, xmax=100, colors='k', label=r"True $K_1^c$ = "+str(round(Kc,6)))

for l in L:

    filename = "data/critical_point_L_"+str(l)+"_K_-0.44.txt"
    data.append(np.loadtxt(filename, unpack=True, delimiter=","))
    
for i in range(len(L)):
    
    plt.plot(data[i][0], data[i][2], label="L = "+str(L[i])+", S = "+str(round(L[i]/b)), c=colors[i], linewidth=2)
    
    n = int(50*(np.log2(L[i])-1))
    avg = np.average(data[i][2][-n:])
    plt.hlines(avg, xmin=-10, xmax=100, linestyle="--", colors=colors[i], label=r"$K_1^c$ = "+str(round(avg,6)))
    
    print("%i %.15f" % (L[i], avg))
    
plt.hlines(Kc, xmin=-10, xmax=100, colors='k', linestyle="--", label=r"Avg $K_1^c$ over last 50 iterations")

plt.title(r"Locating critical nearest neighbor coupling $K_1^c$")
plt.ylabel(r"$K_1^c$")
plt.xlabel(r"Iteration")
plt.xlim(0, 85)
plt.legend(ncol=2)
plt.savefig("critical_point.pdf", format="pdf")
