import numpy as np
import matplotlib.pyplot as plt
import os
import glob

# returns list of training files in data directory
def get_filenames():

    format = 'data/train*.txt'
    filenames = glob.glob(format)
    print('Found %i training files.' % len(filenames))
    
    return filenames

# returns N and T extracted from file name
def get_params(filename):
    
    start = filename.find('_N')+3
    end = filename.find('_', start)
    N = int(filename[start:end])
    
    start = filename.find('_T')+3
    end = filename.find('.txt', start)
    T = float(filename[start:end])
    
    return N, T


# returns final weights extracted from end of file
def get_final_weights(filename):

    file = open(filename, "r")
    weights = []
    
    for line in file:
        if "Final Weights" in line:
            start = line.find(':')+1
            weights_str = line[start:].split()
            for w in weights_str:
                weights.append(float(w))
    
    return weights


Tc = 2/np.log(1+np.sqrt(2))


for filename in get_filenames():
    
    N, T = get_params(filename)
    weights = get_final_weights(filename)
    
    if N == 4:
        c = 'r'
        
    if N == 8:
        c = 'g'
        
    if N == 16:
        c = 'b'
    
    temp = (T/Tc)*np.ones(len(weights))
    plt.scatter(temp, weights, color=c)

plt.show()
