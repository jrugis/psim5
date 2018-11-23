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

## get the time values associated with the saved data
def get_time_vals(fname):
  f = open(fname + ".dat", "r") # get the saved data stride
  for line in f:
    if "Tstride" in line:
      tstride = int(line.rstrip().split()[1])
      break
  f.close()
  tvals = [] # list of time values
  n = int(0)
  f = open(fname + ".out", "r") 
  for line in f: 
    if "<Acinus> t:" in line: # iterate through the saved time step values
      if(float(line.rstrip().split()[2]) == 0.0): continue  # skip start time 0.0
      n += 1
      if (n % tstride) == 0:
        tvals.append(float(line.rstrip().split()[2])) # add to list of time values
  f.close()
  return np.asarray(tvals) # return time values in a numpy array

# get cell node count
def get_node_count(fname):
  f = open(fname, "r")
  for line in f:
    if line.startswith("<CellMesh> vertices_count:"):
      v = line.rstrip().split()
  f.close()
  return int(v[2])


##################################################################
# main program
##################################################################

print("create summary plot")

fig, plots = plt.subplots(len(dtypes), ncells, sharex='col', squeeze=False)
plt.subplots_adjust(wspace = 0.5)
fig.set_size_inches(ncells * 3.8, len(dtypes) * 2.5)
fig.text(0.02, 0.96, os.getcwd(), fontsize=10)

x = get_time_vals("a1") # get the x-axis time values
for cell in range(ncells):
  dname = "a1c" + str(cell + 1)
  nodes = get_node_count(dname + ".out")
  plots[len(dtypes)-1, cell].set_xlabel(" time (s)")
  plots[0, cell].set_title("Cell " + str(cell+1))
  for i, dtype in enumerate(dtypes):
    data = get_data(dname + '_' + dtype + ".bin", nodes, x.shape[0])
    y1 = data[0]
    y2 = data[1]
    plots[i, cell].plot(x, y1, x, y2, color='steelblue')
    plots[i, cell].fill_between(x, y1, y2, color='steelblue')
    if cell == 0:
      plots[i, cell].set_ylabel(dtype + " (uM)")

fig.savefig("summary_plot.pdf")


