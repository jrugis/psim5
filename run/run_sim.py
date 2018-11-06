from __future__ import print_function
import os
import time
import subprocess

##################################################################
# main program
##################################################################

parms = "parms_test_0.dat"
mesh_base = "4sim_out_N4_p3-p2-p4-Xtet.bin"

print("psim5:", parms)

# create the top level results directory
csdir = os.getcwd()
if not os.path.exists("results"):
    os.mkdir("results")
path = "results/" + time.strftime("%y%m%d_%H%M%S")
os.mkdir(path)
os.chdir(path)

# create parms directory
path = parms.split('.')[0]
os.mkdir(path)
os.chdir(path)

# copy some files into results directory
os.system("cp " + csdir + "/psim5 .")
os.system("chmod 770 psim5")
os.system("cp " + csdir + "/run.sl .")
os.system("cp " + csdir + "/summary_plot.py .")
os.system("cp " + csdir + "/" + parms + " a1.dat")
for cell in range(1, 8):
  mesh = mesh_base.replace('X', str(cell))
  
  # first check if the mesh exists (doesn't for all percentage lumen)
  if not os.path.exists(csdir + "/meshes/" + mesh):
    print("Skipping cell {0} as no mesh file exists ({1})!".format(cell, mesh))
    continue

  # copy data files and execute the simulation
  os.system("cp " + csdir + "/meshes/" + mesh + " a1c" + str(cell) + ".bmsh")
  os.system("cp " + csdir + "/" + parms + " a1c" + str(cell) + ".dat")

# submit slurm script
job_output = subprocess.check_output("sbatch run.sl", stderr=subprocess.STDOUT, shell=True)

# go back to top level
os.chdir(csdir)