from __future__ import print_function
import os
import time
import subprocess
import sys

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))

##################################################################
# functions
##################################################################

# create two parameter lists from a sweep parameters file
def get_sweep_parms(fname):
  p1 = []
  p2 = []
  first_parm = ""
  second_parm = ""
  for line in open(fname, "r"):
    if line[0] == '#': continue         # comment line
    if first_parm == "":                # first paramter
      first_parm = line.split()[0]
      p1.append(line)
      continue
    if line.split()[0] == first_parm:
      p1.append(line)
      continue
    if second_parm == "":               # second paramter
      second_parm = line.split()[0]
      p2.append(line)
      continue
    if line.split()[0] == second_parm:
      p2.append(line)
      continue
    print("error: too many sweep parameters")
    quit()
  if len(p1) == 0:
    print("error: no sweep parameters")
    quit()
  if len(p2) == 0: p2.append("")
  return p1, p2

# make a directory-naming-friendly label from a parameter string
def make_label(s):
  tokens = s.strip().replace(".", "p").split()
  return "-".join(tokens)

# replace line(s) in a parameters file
def replace_line(fname, s1, s2):
  s1_found = (s1 == "")
  s2_found = (s2 == "")
  if (s1_found & s2_found): return  # do nothing if no parameters

  s1_parm = s1.strip().split()[0] # first parameter name
  if not s2_found: s2_parm = s2.strip().split()[0] # second parameter name

  ftemp_name = "temp.dat"
  ftemp = open(ftemp_name, "w")
  for line in open(fname, "r"):
    line_parm = line.strip().split()[0]
    if not s1_found:
      if s1_parm == line_parm:
        ftemp.write(s1) # replace line
        s1_found = True
        continue
    if not s2_found:
      if s2_parm == line_parm:
        ftemp.write(s2) # replace line
        s2_found = True
        continue
    ftemp.write(line) # unchanged line
  ftemp.close()
  os.rename(ftemp_name, fname)

  if s1_found == False:           # error checking
    print("error: no such sweep parameter", s1_parm)
    quit()
  if s2_found == False:           # error checking
    print("error: no such sweep parameter", s2_parm)
    quit()
  return

def replace_script_dir(fname):
  with open(fname) as fh:
    lines = fh.readlines()
  for i, line in enumerate(lines):
    if "SCRIPT_DIR" in line:
      lines[i] = line.replace("SCRIPT_DIR", SCRIPT_DIR)
  with open(fname, "w") as fh:
    fh.write("".join(lines))

##################################################################
# main program
##################################################################

print("psim5")
run_dir = os.getcwd()

if(len(sys.argv) < 7):
  print("error: missing argument(s)")
  print("usage: python run_sim.py <slurm-file> <flow-parameter-file> <flow-initial-conditions-file> <flow-adjacency-matrix-file> <calcium-parameter-file> <optional-parameter-sweep-file>")
  quit()

slurm = sys.argv[1] # slurm file
if not os.path.exists(run_dir + "/" + slurm):
  print("no such slurm file: " + slurm)
  quit()

fparms = sys.argv[2] # flow parameters file
if not os.path.exists(run_dir + "/" + fparms):
  print("no such parameter file: " + fparms)
  quit()

finit = sys.argv[3] # flow initial conditions file
if not os.path.exists(run_dir + "/" + fparms):
  print("no such parameter file: " + fparms)
  quit()

fadj = sys.argv[4] # flow adjacency matrix file
if not os.path.exists(run_dir + "/" + fparms):
  print("no such parameter file: " + fparms)
  quit()

parms = sys.argv[5] # calcium parameters file
if not os.path.exists(run_dir + "/" + parms):
  print("no such parameter file: " + parms)
  quit()

p1_array = [""]
p2_array = [""]
if(len(sys.argv) >= 7):
  sweep = sys.argv[6] # parameter sweep file
  if not os.path.exists(run_dir + "/" + sweep):
    print("no such parameter sweep file: " + sweep)
    quit()
  p1_array, p2_array = get_sweep_parms(sweep)

mesh_base = "4sim_out_N4_p3-p2-p4-Xtet.bin"

# create the top level results directory
list = run_dir.split("/")[-4:]
results_dir = "/nesi/nobackup/" + list[0] + "/" + list[1] + "/" + list[2] + "/" + list[3]
results_dir += "/results/" + time.strftime("%y%m%d_%H%M%S")
os.system("mkdir -p " + results_dir)
os.chdir(results_dir)

# setup parameter sweep directories
f1 = open("dirs.txt", "w") # create parameters directory list file
for p1 in p1_array:
  for p2 in p2_array:
    # create parameter directory
    parm_dir = parms.split('.')[0] + '-' + make_label(p1) + '-' + make_label(p2)
    os.mkdir(parm_dir)
    os.chdir(parm_dir)
    f1.write(parm_dir + "\n")

    # copy some files into parameter directory
    os.system("cp " + run_dir + "/psim5 .")
    os.system("chmod 770 psim5")
    os.system("cp " + run_dir + "/" + slurm + " ../run.sl")
    replace_script_dir("../run.sl")
    os.system("cp " + run_dir + "/" + finit + " flow_init.dat")
    os.system("cp " + run_dir + "/" + fadj + " flow_adj.dat")
    os.system("cp " + run_dir + "/" + fparms + " l1.dat")
    os.system("cp " + run_dir + "/" + parms + " a1.dat")
    replace_line("a1.dat", p1, p2)
    for cell in range(1, 8): #copy mesh files
      mesh = mesh_base.replace('X', str(cell))
      if not os.path.exists(run_dir + "/meshes/" + mesh):
        print("no such mesh file: " + mesh)
        quit()
      os.system("cp " + run_dir + "/meshes/" + mesh + " a1c" + str(cell) + ".bmsh")
      # copy matlab dfa
      dfafn = "matlab_dfa_cell%d.dat" % cell
      if not os.path.exists(run_dir + "/" + dfafn):
        print("no such distance file: " + dfafn)
        quit()
      os.system("cp " + run_dir + "/" + dfafn + " .")
      # copy matlab dfb
      dfbfn = "matlab_dfb_cell%d.dat" % cell
      if not os.path.exists(run_dir + "/" + dfbfn):
        print("no such distance file: " + dfafn)
        quit()
      os.system("cp " + run_dir + "/" + dfbfn + " .")
      # copy matlab apical triangles
      apicalfn = "matlab_apical_cell_%d_zeroidx.dat" % cell
      if not os.path.exists(run_dir + "/" + apicalfn):
        print("no such distance file: " + apicalfn)
        quit()
      os.system("cp " + run_dir + "/" + apicalfn + " .")
      # copy matlab basal triangles
      basalfn = "matlab_basal_cell_%d_zeroidx.dat" % cell
      if not os.path.exists(run_dir + "/" + basalfn):
        print("no such distance file: " + basalfn)
        quit()
      os.system("cp " + run_dir + "/" + basalfn + " .")

    os.chdir("..")
f1.close()

# submit slurm script for an array job
cmd = "sbatch --array=1-" + str(len(p1_array) * len(p2_array)) + " run.sl"
job_output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)

# go back to top level
os.chdir(run_dir)
