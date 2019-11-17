#!/bin/bash -e
#SBATCH --job-name=psim5           # job name (shows up in the queue)
#SBATCH --account=nesi00119        # Project Account
#SBATCH --time=01:00:00            # Walltime (HH:MM:SS)
#SBATCH --mem-per-cpu=4000         # memory/cpu (in MB, set to half of what's actually required)
##SBATCH --partition=prepost        # 3 hours, 2 (36) cores,      15GB
#SBATCH --partition=large          # 3 days,  1024 (8424) cores, 3GB
##SBATCH --partition=bigmem         # 7 days,  108 (108) cores,   15GB
#SBATCH --hint=nomultithread       # don't use hyperthreading
#SBATCH --ntasks=8                 # number of tasks (e.g. MPI)
#SBATCH --switches=1@0:10:00       # run on one infiniband switch, wait for 10 minutes

ml SUNDIALS/4.1.0-gimkl-2018b

echo $HOSTNAME
echo "task array id: $SLURM_ARRAY_TASK_ID"

# directory associated with job array
job_dir=$( head -n $SLURM_ARRAY_TASK_ID dirs.txt | tail -1 )
echo $job_dir

cd $job_dir
srun --ntasks=8 psim5
rm psim5

ml Python/3.7.3-gimkl-2018b
srun --ntasks=1 python "SCRIPT_DIR/summary_plot.py"
srun --ntasks=1 python "SCRIPT_DIR/summary_plot_averages.py"
