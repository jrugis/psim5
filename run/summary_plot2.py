from __future__ import print_function

import matplotlib as mpl
mpl.use('Agg') # because there's no X11 display

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import struct
import ctypes
import time
import argparse


# parse command line arguments
parser = argparse.ArgumentParser(description="Create summary plot of psim5 simulation")
parser.add_argument("--no-ca", action="store_true", help="Don't plot ca")
parser.add_argument("--no-ip3", action="store_true", help="Don't plot ip3")
parser.add_argument("--no-volume", action="store_true", help="Don't plot volume")
parser.add_argument("--no-ffr", action="store_true", help="Don't plot fluid flow rate")
parser.add_argument("--ncells", type=int, default=7, help="Number of cells (default is 7)")
parser.add_argument("--nintra", type=int, default=8, help="Number of intracellular variables (default is 8)")
parser.add_argument("--font-size", type=int, default=16, help="Font size for matplotlib (default is 16)")
parser.add_argument("-o", "--output", default="summary_plot2.pdf", help="File name to save plot to (default is summary_plot2.pdf)")
args = parser.parse_args()

plt.rcParams.update({'font.size': args.font_size})
ncells = args.ncells
nintra = args.nintra
dtypes = []
if not args.no_ca:
    dtypes.append("ca")
if not args.no_ip3:
    dtypes.append("ip3")
if not args.no_volume:
    dtypes.append("volume")
plot_ffr = False if args.no_ffr else True

##################################################################
# ctypes lib for loading data
##################################################################

script_dir = os.path.abspath(os.path.dirname(__file__))
lib_path = os.path.join(script_dir, "_summary_plot.so")
if not os.path.exists(lib_path):
    print("you must run 'make'")
    sys.exit(1)
_lib = ctypes.CDLL(lib_path)
_lib.load_summary_plot_data.restype = None
_lib.load_summary_plot_data.argtypes = [
  ctypes.c_char_p,
  ctypes.c_int,
  ctypes.c_int,
  np.ctypeslib.ndpointer(dtype=np.float32),
  np.ctypeslib.ndpointer(dtype=np.float32),
  np.ctypeslib.ndpointer(dtype=np.float32),
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
#fig.set_size_inches(ncells * 3.8, nplots * 2.5)
fig.set_size_inches(ncells * 7.6, nplots * 5.0)
fig.text(0.02, 0.96, os.getcwd(), fontsize=20)

# get the x-axis time values
x = get_time_vals("a1")

# load fluid flow results
if plot_ffr or "volume" in dtypes:
    volumes, ffr = get_data_fluid_flow(x.shape[0])

# main loop over plots
for cell in range(ncells):
  dname = "a1c" + str(cell + 1)
  nodes = get_node_count(dname + ".out")
  if cell == 0 and plot_ffr:
    pass
  else:
    plots[len(dtypes)-1, cell].set_xlabel(" time (s)")
  plots[0, cell].set_title("Cell " + str(cell+1))

  # load ca and ip3 data
  if "ca" in dtypes or "ip3" in dtypes:
      ca_apical = np.empty(x.shape[0], np.float32)
      ca_basal = np.empty(x.shape[0], np.float32)
      ip_apical = np.empty(x.shape[0], np.float32)
      ip_basal = np.empty(x.shape[0], np.float32)
      _lib.load_summary_plot_data(dname.encode("utf-8"), cell + 1, x.shape[0], ca_apical, ca_basal, ip_apical, ip_basal)

  # plot ca
  if "ca" in dtypes:
      ca_index = dtypes.index("ca")
      plots[ca_index, cell].plot(x, ca_apical, color='blue', label='apical')
      plots[ca_index, cell].plot(x, ca_basal, color='red', label='basal')
      plots[ca_index, cell].legend(loc='best')

  # plot ip3
  if "ip3" in dtypes:
      ip_index = dtypes.index("ip3")
      plots[ip_index, cell].plot(x, ip_apical, color='blue', label='apical')
      plots[ip_index, cell].plot(x, ip_basal, color='red', label='basal')
      plots[ip_index, cell].legend(loc='best')

  # plot volumes
  if "volume" in dtypes:
      vol_index = dtypes.index("volume")
      plots[vol_index, cell].plot(x, volumes[:, cell], color='blue')

for i, dtype in enumerate(dtypes):
    if dtype == 'volume':
        plots[i, 0].set_ylabel(dtype + " ($\mu$m$^3$)")
    else:
        plots[i, 0].set_ylabel(dtype + " ($\mu$M)")

if plot_ffr:
  for i in range(1, ncells):
    plots[-1, i].axis('off')
    plots[-2, i].xaxis.set_tick_params(which='both', labelbottom=True)
  plots[-1, 0].set_title("Total fluid flow rate")
  plots[-1, 0].plot(x, ffr, color='blue')
  plots[-1, 0].set_xlabel(" time (s)")
  plots[-1, 0].set_ylabel(" Flow rate ($\mu m^3$/sec)")

fig.savefig(args.output)
