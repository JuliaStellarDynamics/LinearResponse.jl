"""
InspectLinearResponse.py

a simple Python script to read output files from CallAResponse

"""

# standard python modules
import numpy as np

# you will need an HDF5 reader
import h5py

# read in file structures
import glob

# plotting utilities
import matplotlib.pyplot as plt;import matplotlib as mpl;import matplotlib.cm as cm;import matplotlib.colors as colors;cmap = mpl.cm.inferno

# set some basic parameters for plots
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 0.75
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 0.75
mpl.rcParams['ytick.minor.visible'] = True

# set the model parameters: draw these from your input file
model = "PlummerE"
dfname = "roiinf"
harmonic = 1
n1max = 5
rb = 20.0
Ku = 202
Kv = 200

# set directory structure
wmatdir  = "wmat/"
gfuncdir = "gfunc/"
modedir  = "xifunc/"

# let's skip straight to the results: the mode shape
filename = modedir+"ModeShape_{}_df_{}_l_{}_n1_{}_rb_{}_Ku_{}_Kv_{}.h5".format(model,dfname,harmonic,n1max,rb,Ku,Kv)
f = h5py.File(filename, 'r')
print(f.keys())

plt.figure(figsize=(4.5,3),facecolor='white')
modeval = (f['ModePotentialShape'][:])
plt.plot(f['ModeRadius'][:],modeval/np.nanmax(modeval),color='black')
plt.xlabel('Radius')
plt.ylabel('Amplitude')
plt.title("ModeShape_{}_df_{}_l_{}_n1_{}_rb_{}_Ku_{}_Kv_{}.h5".format(model,dfname,harmonic,n1max,rb,Ku,Kv),size=8)
plt.tight_layout()
plt.savefig('figures/ModeShape.png')
