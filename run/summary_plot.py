from __future__ import print_function

import matplotlib as mpl
mpl.use('Agg') # because there's no X11 display

import matplotlib.pyplot as plt
import numpy as np
import os
import struct

ncells = 7
dtypes = ["ca", "ip3", "cer"]

##################################################################
# functions
##################################################################

## read in a simulation data file
def get_data(fname, rows, cols):
  f1 = open(fname, "rb")
  data = np.zeros((2, cols), dtype=np.float32)
  buf = np.zeros(rows, dtype=np.float32)
  for c in range(cols):     # the data is in column order
    for r in range(rows):
      buf[r] = struct.unpack('f', f1.read(4))[0]
    data[0, c] = buf.min()
    data[1, c] = buf.max()
  f1.close()
  return data

# get simulation time data
# NOTE: adjusted for possible time step stride
def get_sim_time(fname):
  f1 = open(fname, "r")
  for line in f1:
    if line.startswith("%! delT totalT Tstride"):
      v = next(f1).rstrip().split()
  f1.close()
  delt = float(v[0])
  totalt = float(v[1])
  tstride = float(v[2])
  steps = int(totalt / delt) 
  ndelt = delt * tstride
  nsteps = int(steps / tstride)
  ntotalt = ndelt * nsteps
  return ndelt, ntotalt, nsteps

# get cell node count
def get_node_count(fname):
  f1 = open(fname, "r")
  for line in f1:
    if line.startswith("<CellMesh> vertices_count:"):
      v = line.rstrip().split()
  f1.close()
  return int(v[2])


##################################################################
# main program
##################################################################

print("create summary plot")

fig, plots = plt.subplots(len(dtypes), ncells, sharex='col', squeeze=False)
plt.subplots_adjust(wspace = 0.5)
fig.set_size_inches(ncells * 3.8, len(dtypes) * 2.5)
fig.text(0.02, 0.96, os.getcwd(), fontsize=10)

delt, totalt, tsteps = get_sim_time("a1.dat")
for cell in range(ncells):
  dname = "a1c" + str(cell + 1)
  nodes = get_node_count(dname + ".out")
  plots[len(dtypes)-1, cell].set_xlabel(" time (s)")
  plots[0, cell].set_title("Cell " + str(cell+1))
  for i, dtype in enumerate(dtypes):
    data = get_data(dname + '_' + dtype + ".bin", nodes, tsteps)
    x = np.linspace(delt, totalt, tsteps)
    y1 = data[0]
    y2 = data[1]
    plots[i, cell].plot(x, y1, x, y2, color='steelblue')
    plots[i, cell].fill_between(x, y1, y2, color='steelblue')
    if cell == 0:
      plots[i, cell].set_ylabel(dtype + " (uM)")

fig.savefig("summary_plot.pdf")


