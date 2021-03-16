import numpy as np
import matplotlib.pyplot as plt
from tools import get_params, get_results, get_filenames


Tc = 2/np.log(1+np.sqrt(2))

for filename in get_filenames():
    
    N, T = get_params(filename)
    results = get_final_weights(filename, 10000)
    
    if N == 16:
        c = 'b'
        
    plt.scatter(T/Tc, results[0]/Tc, c=c)
    
    
temp = np.linspace(0.5, 1.5)
plt.plot(temp, temp, c='k')

plt.xlabel(r"$T/T_c$")
plt.ylabel(r"$T_{pred}/T_c$")
plt.grid(alpha=0.2)
    

plt.show()

