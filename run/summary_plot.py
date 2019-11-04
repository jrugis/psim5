from __future__ import print_function

import matplotlib as mpl
mpl.use('Agg') # because there's no X11 display

import matplotlib.pyplot as plt
import numpy as np
import os
import struct
import ctypes
import time

ncells = 7
nintra = 8
dtypes = ["ca", "ip3", "cer", "volume"]
plot_ffr = True


##################################################################
# ctypes lib for loading data
##################################################################

script_dir = os.path.abspath(os.path.dirname(__file__))
try:
    _lib = ctypes.CDLL(os.path.join(script_dir, "_get_data.so"))
except OSError:
    _lib = None
    print("Using python get_data - compile C library for better performance")
else:
    _lib.get_data.restype = None
    _lib.get_data.argtypes = [
      ctypes.c_char_p,
      ctypes.c_int,
      ctypes.c_int,
      np.ctypeslib.ndpointer(dtype=np.float32),
    ]

##################################################################
# functions
##################################################################

def get_data_fluid_flow(ntime):
  fname = "l1_results.dat"
  data = np.loadtxt(fname)
  assert data.shape[0] == ntime, "Wrong number of lines in fluid flow results"

  vol_cols = [i * nintra + 0 for i in range(ncells)]
  volumes = data[:, vol_cols]
  ffr = data[:, -1]

  return volumes, ffr

def get_data(fname, rows, cols):
  if _lib is None:
    data = get_data_py(dname + '_' + dtype + ".bin", nodes, x.shape[0])
  else:
    data = get_data_c(dname + '_' + dtype + ".bin", nodes, x.shape[0])
  return data

## read in a simulation data file
def get_data_c(fname, rows, cols):
  data = np.zeros((2, cols), dtype=np.float32)
  _lib.get_data(fname.encode("utf-8"), rows, cols, data)
  return data

## read in a simulation data file
def get_data_py(fname, rows, cols):
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

nplots = len(dtypes) + 1 if plot_ffr else len(dtypes)
fig, plots = plt.subplots(nplots, ncells, sharex='col', squeeze=False)
plt.subplots_adjust(wspace = 0.5)
fig.set_size_inches(ncells * 3.8, nplots * 2.5)
fig.text(0.02, 0.96, os.getcwd(), fontsize=10)

x = get_time_vals("a1") # get the x-axis time values
volumes, ffr = get_data_fluid_flow(x.shape[0])
for cell in range(ncells):
  dname = "a1c" + str(cell + 1)
  nodes = get_node_count(dname + ".out")
  if cell == 0 and plot_ffr:
    pass
  else:
    plots[len(dtypes)-1, cell].set_xlabel(" time (s)")
  plots[0, cell].set_title("Cell " + str(cell+1))
  for i, dtype in enumerate(dtypes):
    if dtype == 'volume':
      plots[i, cell].plot(x, volumes[:, cell], color='steelblue')
      if cell == 0:
        plots[i, cell].set_ylabel(dtype + " ($\mu m^3$)")
    else:
      data = get_data(dname + '_' + dtype + ".bin", nodes, x.shape[0])
      y1 = data[0]
      y2 = data[1]
      plots[i, cell].plot(x, y1, x, y2, color='steelblue')
      plots[i, cell].fill_between(x, y1, y2, color='steelblue')
      if cell == 0:
        plots[i, cell].set_ylabel(dtype + " (uM)")

if plot_ffr:
  for i in range(1, ncells):
    plots[-1, i].axis('off')
  plots[-1, 0].set_title("Total fluid flow rate")
  plots[-1, 0].plot(x, ffr, color='steelblue')
  plots[-1, 0].set_xlabel(" time (s)")
  plots[-1, 0].set_ylabel(" Flow rate ($\mu m^3$/sec)")

fig.savefig("summary_plot.pdf")
