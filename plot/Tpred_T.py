import numpy as np
import matplotlib.pyplot as plt
from tools import get_params, get_results, get_filenames


Tc = 2/np.log(1+np.sqrt(2))

for filename in get_filenames():
    
    N, T = get_params(filename)
    results = get_results(filename, 10000)
    
    if N == 8:
        c = 'r'
    
    if N == 16:
        c = 'b'
    
        
    plt.scatter(T/Tc, results[0]/Tc, c=c, alpha=0.7, edgecolor='k', linewidth=0.5)
    
plt.scatter(0, 0, c='r', label='N = '+str(8), alpha=0.7, edgecolor='k', linewidth=0.5)
plt.scatter(0, 0, c='b', label='N = '+str(16), alpha=0.7, edgecolor='k', linewidth=0.5)
temp = np.linspace(0.5, 1.5)
plt.plot(temp, temp, c='k', linestyle='dashed')

plt.xlim(0.49, 1.51)
plt.xlabel(r"$T/T_c$")
plt.ylabel(r"$T_{pred}/T_c$")
plt.grid(alpha=0.2)
plt.legend()

plt.savefig("plot/Tpred_T.pdf", format='pdf')

