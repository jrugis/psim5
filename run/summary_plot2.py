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
parser.add_argument("--cells", type=int, nargs='*', default=[1, 2, 3, 4, 5, 6, 7], help="Cells to plot (default is 1 2 3 4 5 6 7)")
parser.add_argument("--nintra", type=int, default=8, help="Number of intracellular variables (default is 8)")
parser.add_argument("--font-size", type=int, default=16, help="Font size for matplotlib (default is 16)")
parser.add_argument("-o", "--output", default="summary_plot2.pdf", help="File name to save plot to (default is summary_plot2.pdf)")
args = parser.parse_args()

plt.rcParams.update({'font.size': args.font_size})
ncells = len(args.cells)
print("plotting cells:", args.cells)
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
_lib.load_summary_plot_data.restype = ctypes.c_int
_lib.load_summary_plot_data.argtypes = [
  ctypes.c_char_p,
  ctypes.c_int,
  ctypes.c_int,
  np.ctypeslib.ndpointer(dtype=np.float32),
  np.ctypeslib.ndpointer(dtype=np.float32),
  ctypes.c_int,
  np.ctypeslib.ndpointer(dtype=np.int32),
  ctypes.c_int,
  np.ctypeslib.ndpointer(dtype=np.int32),
]

##################################################################
# functions
##################################################################

def get_data_fluid_flow(ntime):
  fname = "l1_results.dat"
  if not os.path.exists(fname):
      return None, None

  data = np.loadtxt(fname)
  assert data.shape[0] == ntime, "Wrong number of lines in fluid flow results"

  vol_cols = [i * nintra + 0 for i in args.cells]
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

# get the x-axis time values
x = get_time_vals("a1")

# load fluid flow results
if plot_ffr or "volume" in dtypes:
    volumes, ffr = get_data_fluid_flow(x.shape[0])
    # if fluid flow was disabled don't plot it
    if volumes is None and ffr is None:
        print("Warning: disabling fluid flow and volume plotting since no lumen results file exists")
        plot_ffr = False
        if "volume" in dtypes:
            dtypes.remove("volume")

print("create summary plot")

nplots = len(dtypes) + 1 if plot_ffr else len(dtypes)
fig, plots = plt.subplots(nplots, ncells, sharex='col', squeeze=False)
plt.subplots_adjust(wspace = 0.5)
#fig.set_size_inches(ncells * 3.8, nplots * 2.5)
fig.set_size_inches(ncells * 7.6, nplots * 5.0)
fig.text(0.02, 0.96, os.getcwd(), fontsize=20)

# main loop over plots
for celli, cell in enumerate(args.cells):
  dname = "a1c" + str(cell)
  print("plotting {}".format(dname))
  nodes = get_node_count(dname + ".out")
  if celli == 0 and plot_ffr:
    pass
  else:
    plots[len(dtypes)-1, celli].set_xlabel(" time (s)")
  plots[0, celli].set_title("Cell " + str(cell))

  # load ca and ip3 data
  if "ca" in dtypes or "ip3" in dtypes:
      # first we load the list of apical and basal nodes
      apical_nodes = np.fromfile("apical_nodes_{}.bin".format(dname), dtype=np.int32)
      basal_nodes = np.fromfile("basal_nodes_{}.bin".format(dname), dtype=np.int32)
      print("  num apical nodes = {0}, num basal nodes = {1}".format(len(apical_nodes), len(basal_nodes)))

      # now load the result and calculate averages at apical and basal nodes
      if "ca" in dtypes:
          # for calcium
          ca_apical = np.empty(x.shape[0], np.float32)
          ca_basal = np.empty(x.shape[0], np.float32)
          filename = dname + "_ca.bin"
          status = _lib.load_summary_plot_data(filename.encode("utf-8"), x.shape[0], nodes, ca_apical, ca_basal,
                  len(apical_nodes), apical_nodes, len(basal_nodes), basal_nodes)
          if status < 0:
              sys.exit("Error open result file: {}".format(result_file))
          if status > 0:
              sys.exit("Error reading row in result file: {}".format(status))

      if "ip3" in dtypes:
          # for ip3
          ip_apical = np.empty(x.shape[0], np.float32)
          ip_basal = np.empty(x.shape[0], np.float32)
          filename = dname + "_ip3.bin"
          status = _lib.load_summary_plot_data(filename.encode("utf-8"), x.shape[0], nodes, ip_apical, ip_basal,
                  len(apical_nodes), apical_nodes, len(basal_nodes), basal_nodes)
          if status < 0:
              sys.exit("Error open result file: {}".format(result_file))
          if status > 0:
              sys.exit("Error reading row in result file: {}".format(status))

  # plot ca
  if "ca" in dtypes:
      ca_index = dtypes.index("ca")
      plots[ca_index, celli].plot(x, ca_apical, color='blue', label='apical')
      plots[ca_index, celli].plot(x, ca_basal, color='red', label='basal')
      plots[ca_index, celli].legend(loc='best')

  # plot ip3
  if "ip3" in dtypes:
      ip_index = dtypes.index("ip3")
      plots[ip_index, celli].plot(x, ip_apical, color='blue', label='apical')
      plots[ip_index, celli].plot(x, ip_basal, color='red', label='basal')
      plots[ip_index, celli].legend(loc='best')

  # plot volumes
  if "volume" in dtypes:
      vol_index = dtypes.index("volume")
      plots[vol_index, celli].plot(x, volumes[:, celli], color='blue')

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
