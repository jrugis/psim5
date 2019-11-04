from __future__ import print_function

import os

import matplotlib as mpl
mpl.use('Agg') # because there's no X11 display

import matplotlib.pyplot as plt
import numpy as np


NUMCELLS = 7
NUMINTRAVARS = 8


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


##################################################################
# main program
##################################################################

print("create volumes plot")

#fig, plots = plt.subplots(len(dtypes), ncells, sharex='col', squeeze=False)
#plt.subplots_adjust(wspace = 0.5)
#fig.set_size_inches(ncells * 3.8, len(dtypes) * 2.5)
#fig.text(0.02, 0.96, os.getcwd(), fontsize=10)

fig = plt.figure()
dirtext = os.path.sep.join(os.getcwd().split(os.path.sep)[-2:])
fig.text(0.02, 0.96, dirtext, fontsize=10)
ax = fig.add_subplot(1, 1, 1)

# get the x-axis time values
x = get_time_vals("a1")

# volume columns (volume is first intracellular variable)
vol_cols = [i * NUMINTRAVARS + 0 for i in range(NUMCELLS)]
ion_res = np.loadtxt("l1_results.dat", skiprows=0)
volumes = ion_res[:, vol_cols]

for cell in range(NUMCELLS):
    y = volumes[:, cell]
    ax.plot(x, y, label="Cell {}".format(cell + 1))

ax.legend(loc='best')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Volume ($\mu$m$^3$)')

fig.savefig("volumes_plot.pdf")
