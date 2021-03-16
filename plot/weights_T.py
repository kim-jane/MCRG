import numpy as np
import matplotlib.pyplot as plt
from tools import get_params, get_final_weights, get_filenames


Tc = 2/np.log(1+np.sqrt(2))
#colors = ['indianred', 'orange', 'forestgreen', 'dodgerblue']

for filename in get_filenames():
    
    N, T = get_params(filename)
    weights = get_final_weights(filename)
    
    if N == 8:
        c = 'r'
    
    if N == 16:
        c = 'b'
    
    plt.scatter( (T/Tc)*np.ones(len(weights)), weights, c=c, alpha=0.7, edgecolor='k', linewidth=0.5)
    

plt.scatter(0, 0, c='r', label='N = '+str(8), alpha=0.7, edgecolor='k', linewidth=0.5)
plt.scatter(0, 0, c='b', label='N = '+str(16), alpha=0.7, edgecolor='k', linewidth=0.5)

plt.ylim(-1,1)
plt.xlim(0.49, 1.51)
plt.xlabel(r"$T/T_c$")
plt.ylabel(r"Weights")
plt.grid(alpha=0.2)
plt.legend()

plt.savefig("plot/weights_T.pdf", format='pdf')
