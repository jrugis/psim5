from __future__ import print_function
import os
import time
import subprocess
import sys

##################################################################
# main program
##################################################################

if(len(sys.argv) < 3):
  print("error: missing argument(s)")
  print("usage: python run_sim.py <slurm-file-name> <parameter-file-name>")
  quit()

run_dir = os.getcwd()
slurm = sys.argv[1]
if not os.path.exists(run_dir + "/" + slurm):
  print("no such slurm file: " + slurm)
  quit()
parms = sys.argv[2]
if not os.path.exists(run_dir + "/" + parms):
  print("no such parameter file: " + parms)
  quit()
    
print("psim5:", parms)
mesh_base = "4sim_out_N4_p3-p2-p4-Xtet.bin"

# create the top level results directory
list = run_dir.split("/")[-4:]
temp_dir = "/nesi/nobackup/" + list[0] + "/" + list[1] + "/" + list[2] + "/" + list[3] 
temp_dir += "/results/" + time.strftime("%y%m%d_%H%M%S")
os.system("mkdir -p " + temp_dir)
os.chdir(temp_dir)

# create parms directory
path = parms.split('.')[0]
os.mkdir(path)
os.chdir(path)

# copy some files into results directory
os.system("cp " + run_dir + "/psim5 .")
os.system("chmod 770 psim5")
os.system("cp " + run_dir + "/" + slurm + " ./run.sl")
os.system("cp " + run_dir + "/summary_plot.py .")
os.system("cp " + run_dir + "/" + parms + " a1.dat")

for cell in range(1, 8):
  mesh = mesh_base.replace('X', str(cell))
  
  # first check if the mesh exists (doesn't for all percentage lumen)
  if not os.path.exists(run_dir + "/meshes/" + mesh):
    print("Skipping cell {0} as no mesh file exists ({1})!".format(cell, mesh))
    continue

  # copy data files and execute the simulation
  os.system("cp " + run_dir + "/meshes/" + mesh + " a1c" + str(cell) + ".bmsh")

# submit slurm script
job_output = subprocess.check_output("sbatch run.sl", stderr=subprocess.STDOUT, shell=True)

# go back to top level
os.chdir(run_dir)