
# standard python modules
import numpy as np

# plotting utilities
import matplotlib.pyplot as plt;import matplotlib as mpl;import matplotlib.cm as cm;import matplotlib.colors as colors

cmap = mpl.cm.inferno
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 0.75
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 0.75
mpl.rcParams['ytick.minor.visible'] = True

# HDF5 reader for Python
import h5py


# try a wmat view
plt.figure()
f = h5py.File('wmat/wmat_l_2_n1_-2_n2_-2.h5', 'r')

plt.contourf(f['wmat'][:,:,3],cmap=cm.coolwarm)
#plt.contourf(np.log10(f['amat'][:,:]),64,cmap=cm.coolwarm)
#plt.contourf((f['emat'][:,:]),64,cmap=cm.coolwarm)

plt.colorbar()

"""
plt.figure()
f = h5py.File("gfunc/Gfunc_n1_-1_n2_2.200.h5", 'r')
for key in f.keys():
    plt.plot(f[key],color="black")
"""
