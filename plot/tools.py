import numpy as np
import os
import glob


def sort_key(s):
    if '=' in s[0]:
        n = s[0].split('=')[-1]
        return int(n)

# returns list of training files in data directory
def get_filenames(keyword):

    format = 'data/'+keyword+'*.txt'
    filenames = glob.glob(format)
    print('Found %i %s files.' % (len(filenames), keyword))
    
    return filenames

# returns N and T extracted from file name
def get_params(filename):

    start = filename.find('_b')+2
    end = filename.find('_', start)
    b = int(filename[start:end])
    
    start = filename.find('_L')+2
    end = filename.find('_', start)
    L = int(filename[start:end])
    
    start = filename.find('_K')+2
    end = filename.find('.txt', start)
    K = float(filename[start:end])
    
    return b, L, K


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

# returns results at specified cycle number
def get_results(filename, cycle):

    file = open(filename, "r")
    results = []
    
    for line in file:
        
        results_str = line.split()
        
        if len(results_str) > 0:
            if str(cycle) == results_str[0]:
                for r in results_str[1:]:
                    results.append(float(r))
                
    return results


