from __future__ import print_function
import os
import time
import subprocess
import sys

##################################################################
# main program
##################################################################

print("psim5")
run_dir = os.getcwd()

if(len(sys.argv) < 3):
  print("error: missing argument(s)")
  print("usage: python run_sim.py <slurm-file> <parameter-file> <optional-parameter-sweep-file>")
  quit()

slurm = sys.argv[1] # slurm file
if not os.path.exists(run_dir + "/" + slurm):
  print("no such slurm file: " + slurm)
  quit()

parms = sys.argv[2] # parameters file
if not os.path.exists(run_dir + "/" + parms):
  print("no such parameter file: " + parms)
  quit()

if(len(sys.argv) == 4):
  sweep = sys.argv[3] # parameter sweep file
  if not os.path.exists(run_dir + "/" + sweep):
    print("no such parameter sweep file: " + sweep)
    quit()

mesh_base = "4sim_out_N4_p3-p2-p4-Xtet.bin"

# create the top level results directory
list = run_dir.split("/")[-4:]
results_dir = "/nesi/nobackup/" + list[0] + "/" + list[1] + "/" + list[2] + "/" + list[3] 
results_dir += "/results/" + time.strftime("%y%m%d_%H%M%S")
os.system("mkdir -p " + results_dir)
os.chdir(results_dir)

# setup parameter sweep (if any)
f1 = open("temp_dirs.txt", "w") # create parameters directory list file
rows = 1
cols = 1
for r in range(rows):
  for c in range(cols):
    # create parameter directory
    parm_dir = parms.split('.')[0]
    os.mkdir(parm_dir)
    os.chdir(parm_dir)
    f1.write(parm_dir + "\n")

    # copy some files into parameter directory
    os.system("cp " + run_dir + "/psim5 .")
    os.system("chmod 770 psim5")
    os.system("cp " + run_dir + "/" + slurm + " ../run.sl")
    os.system("cp " + run_dir + "/summary_plot.py .")
    os.system("cp " + run_dir + "/" + parms + " a1.dat")
    for cell in range(1, 8): #copy mesh files
      mesh = mesh_base.replace('X', str(cell))
      if not os.path.exists(run_dir + "/meshes/" + mesh):
        print("no such mesh file: " + mesh)
        quit()
      os.system("cp " + run_dir + "/meshes/" + mesh + " a1c" + str(cell) + ".bmsh")
    os.chdir("..")
f1.close()

# submit slurm script
cmd = "sbatch --array=1-" + str(rows * cols) + " run.sl"
job_output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

# go back to top level
os.chdir(run_dir)
