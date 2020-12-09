import numpy as np
import matplotlib.pyplot as plt
import glob

K1 = []
m = []

filelist = glob.glob('data/equilibrate_N_16_*.txt')

for file in filelist:

    i = file.find("K1")+3
    j = file.find("K2")-1
    K1.append(float(file[i:j]))

    data = np.loadtxt(file, delimiter=",")
    m.append(data[-1, 5])
    
    plt.plot(data[:,0], data[:,5])
    
    
plt.legend()
plt.show()

K1, m = zip(*sorted(zip(K1, m)))

plt.figure(figsize=(10,6))
plt.plot(K1, m)
plt.legend()
plt.show()
