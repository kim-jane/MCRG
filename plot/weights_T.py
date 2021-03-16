import numpy as np
import matplotlib.pyplot as plt
from tools import get_params, get_final_weights, get_filenames


Tc = 2/np.log(1+np.sqrt(2))

for filename in get_filenames():
    
    N, T = get_params(filename)
    weights = get_final_weights(filename)
    temp = (T/Tc)*np.ones(len(weights))
    
    if N == 16:
        c = 'b'
        
    plt.scatter(temp, weights, c=c)
    
    
plt.ylim(-1,1)
plt.xlabel(r"$T/T_c$")
plt.ylabel(r"Weights")
plt.grid(alpha=0.2)
    

plt.show()
