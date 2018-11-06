from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import os
import struct
import subprocess
import sys


results_dir = str(sys.argv[1]) # path to results e.g. 181106_011232/parms_test_0
nodes_skip = str(sys.argv[2])  # subsample data stride
steps_skip = str(sys.argv[3])  #

##################################################################
## read in a simulation data file
def get_data(fname, rows, cols):
  f1 = open(fname, "rb")
  data = np.zeros((rows, cols), dtype=np.float32)
  for c in range(cols):     # the data is in column order
    for r in range(rows):
      data[r, c] = struct.unpack('f', f1.read(4))[0]
  f1.close()
  return data

##################################################################
# main program
##################################################################

# create the top level results directory
csdir = os.getcwd()
path = "results/" + results_dir
if not os.path.exists(path):
  os.makedirs(path)
os.chdir(path)

# plot the calcium concentration
path = csdir + "/../run/results/" + results_dir
data = get_data(path + "/a1c2_ca.bin", 13718, 100)[0::int(nodes_skip),0::int(steps_skip)]
p1 = plt.plot(np.transpose(data))
plt.savefig("a1c2_ca.pdf")
plt.show(p1)




# go back to top level
os.chdir(csdir)