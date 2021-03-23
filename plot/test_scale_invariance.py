import numpy as np
import matplotlib.pyplot as plt
from tools import get_params, get_filenames, sort_key

Tc = 2/np.log(1+np.sqrt(2))
lw = 3
c = ['indianred', 'orange', 'forestgreen', 'royalblue']

fig, ax = plt.subplots(figsize=(10,8))

for filename in get_filenames('test_scalar'):

    b, L, K = get_params(filename)
    results = np.loadtxt(filename, unpack=True)
    
    temp = results[1]
    avg_outputL = results[2]
    var_outputL = results[3]
    
    i = int(np.log(L)/np.log(b)-3)
    ax.plot(temp/Tc, avg_outputL, label="L="+str(L), color=c[i], linewidth=lw)
    ax.fill_between(temp/Tc, avg_outputL+np.sqrt(var_outputL), avg_outputL-np.sqrt(var_outputL), color=c[i], alpha=0.2)
    

# sort legend
handles, labels = ax.get_legend_handles_labels()
labels, handles = zip(*sorted(zip(labels, handles), key=sort_key))
ax.legend(handles, labels)

# fig params
plt.ylabel(r'Scalar Output of RGNN')
plt.xlim(temp[0]/Tc, temp[-1]/Tc)
plt.xlabel(r'$T/T_c$')
plt.title('Testing Scale Invariance of RGNN Output')
plt.grid(alpha=0.2)
plt.savefig('test_scale_invariance.pdf', format='pdf')
    
