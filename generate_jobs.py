import numpy as np
import pandas as pd
import itertools
import datetime

dir = "/mnt/home/kimjane7/MCRG"


n_processes = 100
nodes = int(np.ceil(n_processes/20.0))

b = 2
N = [8]

Tc = 2.0/np.log(1.0+np.sqrt(2.0))
T = np.linspace(0.5*Tc, 1.5*Tc, num=5)


for n in N:
    for t in T:
    
        jobfile = "jobs/train_b%d_N%d_T%.7f_np%d.sb"\
                  % (b, n, t, n_processes)
                  
        command = "srun -n %d train %d %d %.15f" \
                  % (n_processes, b, n, t)

        f = open(jobfile, "w")
        f.write("#!/bin/bash --login\n")
        f.write("#SBATCH --time=03:59:00\n")
        f.write("#SBATCH --ntasks=%d\n" % n_processes)
        f.write("#SBATCH --nodes=%d\n" % nodes)
        f.write("#SBATCH --mem-per-cpu=%dG\n" % int(np.ceil(n/8.0)))
        f.write("#SBATCH --job-name='RGNN'\n\n")
        
        f.write("module purge\n")
        f.write("module load GCC/6.4.0-2.28 OpenMPI\n")
        f.write("cd %s\n\n" % dir)
        f.write(command)
        f.close()
