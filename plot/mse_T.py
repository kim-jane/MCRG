import numpy as np
import matplotlib.pyplot as plt
from tools import get_params, get_results, get_filenames


Tc = 2/np.log(1+np.sqrt(2))

fig = plt.figure()
ax = plt.gca()

for filename in get_filenames():
    
    N, T = get_params(filename)
    results = get_results(filename, 10000)
    
    if N == 8:
        c = 'r'
        
    if N == 16:
        c = 'b'
        
    ax.scatter(T/Tc, results[2], c=c, alpha=0.7, edgecolor='k', linewidth=0.5)

plt.scatter(0, 0, c='r', label='N = '+str(8), alpha=0.7, edgecolor='k', linewidth=0.5)
plt.scatter(0, 0, c='b', label='N = '+str(16), alpha=0.7, edgecolor='k', linewidth=0.5)
plt.xlim(0.49, 1.51)
ax.set_yscale('log')
plt.xlabel(r"$T/T_c$")
plt.ylabel(r"MSE")
plt.grid(alpha=0.2)
plt.legend()

plt.savefig("plot/MSE_T.pdf", format='pdf')


