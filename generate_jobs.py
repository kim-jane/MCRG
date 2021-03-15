import numpy as np
import pandas as pd
import itertools
import datetime

dir = "/mnt/home/kimjane7/MCRG"


n_processes = 100
N = [2**i for i in range(3,8)]

Tc = 2.0/np.log(1.0+np.sqrt(2.0))
DeltaT = 0.1
T = np.linspace(Tc-DeltaT, Tc+DeltaT, num=20)

for n in N:
    for t in T:
    
        jobfile = "jobs/train_N%d_T%.7f_np%d.sb"\
                  % (n, t, n_processes)
                  
        command = "srun -n %d train %d %.15f" \
                  % (n_processes, n, t)

        f = open(jobfile, "w")
        f.write("#!/bin/bash --login\n")
        f.write("#SBATCH --time=03:59:00\n")
        f.write("#SBATCH --ntasks=%d\n" % n_processes)
        f.write("#SBATCH --mem-per-cpu=%dG\n" % int(np.ceil(n/8.0)))
        f.write("#SBATCH --job-name='rgnn'\n\n")
        
        f.write("module purge\n")
        f.write("module load GCC/6.4.0-2.28 OpenMPI\n")
        f.write("cd %s\n\n" % dir)
        f.write(command)
        f.close()
