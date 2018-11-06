from __future__ import print_function

import matplotlib as mpl
mpl.use('Agg') # because there's no X11 display

import matplotlib.pyplot as plt
import numpy as np
import os
import struct

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
def get_sim_time(fname):
  f1 = open(fname, "r")
  for line in f1:
    if line.startswith("%! delT totalT"):
      v = next(f1).rstrip().split()
  f1.close()
  return float(v[0]), float(v[1])

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

dname = "a1c2"

delt, totalt = get_sim_time("a1.dat")
tsteps = int(totalt / delt)
nodes = get_node_count(dname + ".out")

# plot the calcium concentration
fig = plt.figure()
data = get_data(dname + "_cer.bin", nodes, tsteps)

x = np.linspace(0.0, totalt, tsteps)
y1 = data[0]
y2 = data[1]
plt.plot(x, y1, x, y2, color='steelblue')
plt.fill_between(x, y1, y2, color='steelblue')
fig.savefig(dname + "_cer.pdf")


