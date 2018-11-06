#!/bin/bash
#SBATCH --job-name=psim5           # job name (shows up in the queue)
#SBATCH --account=nesi00119        # Project Account
#SBATCH --time=01:00:00            # Walltime (HH:MM:SS)
#SBATCH --mem-per-cpu=2000         # memory/cpu (in MB) Should be half?
##SBATCH --partition=prepost        # 3 hours, 4 cores,    15GB
##SBATCH --partition=large          # 3 days,  1024 cores, 3GB
#SBATCH --partition=bigmem         # 3 days,  72 cores,   15GB
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --ntasks=8                 # number of tasks (e.g. MPI)
##SBATCH --ntasks-per-node=8
##SBATCH --ntasks-per-core=8
#SBATCH --switches=1@0:10:00       # run on one infiniband switch, wait for 10 minutes

echo "$HOSTNAME"
srun --ntasks=8 psim5
rm psim5

ml Python/2.7.14-gimkl-2017a
srun --ntasks=1 python "summary_plot.py"
rm summary_plot.py
