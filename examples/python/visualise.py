
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
f1 = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/examples/IsochroneA/wmat/Wmat_IsochroneA_l_2_n1_-1_n2_-2_rb_15.0_Ku_202_Kv_200_Kw_200.h5', 'r')

f2 = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/examples/IsochroneE/wmat/Wmat_IsochroneE_l_2_n1_-1_n2_-2_rb_15.0_Ku_202_Kv_200_Kw_200.h5', 'r')

#>>> f.keys()
#<KeysViewHDF5 ['AEmat', 'BasisParameters', 'ELmat', 'LinearParameters', 'Omgmat', 'OrbitalParameters', 'UVmat', 'jELABmat', 'omgmax', 'omgmin', 'tabvminmax', 'wmat']>

#>>> f1['tabvminmax'].shape
#(202, 2) = (Ku, 2)

plt.contourf(f1['AEmat'][:,:,0]-f2['AEmat'][:,:,0],cmap=cm.coolwarm)


plt.contourf(f1['wmat'][:,:,0],cmap=cm.coolwarm)

plt.contourf(f2['wmat'][:,:,0],cmap=cm.coolwarm)

indx = 2
plt.contourf(f2['wmat'][:,:,indx]-f1['wmat'][:,:,indx],cmap=cm.coolwarm)

#plt.contourf(np.log10(f['amat'][:,:]),64,cmap=cm.coolwarm)
#plt.contourf((f['emat'][:,:]),64,cmap=cm.coolwarm)

# 0 is a, 1 is e
plt.contourf(f2['AEmat'][:,:,1],cmap=cm.coolwarm)

plt.contourf(f1['AEmat'][:,:,1]-f2['AEmat'][:,:,1],cmap=cm.coolwarm)
plt.contourf(f1['AEmat'][:,:,0]-f2['AEmat'][:,:,0],cmap=cm.coolwarm)

# mean is 1.e-8
plt.contourf(f1['Omgmat'][:,:,1]-f2['Omgmat'][:,:,1],cmap=cm.coolwarm)


plt.contourf(f1['ELmat'][:,:,1]-f2['ELmat'][:,:,1],cmap=cm.coolwarm)

plt.contourf(f1['UVmat'][:,:,0]-f2['UVmat'][:,:,0],cmap=cm.coolwarm)
# v accuracy is the same as frequency accuracy?
# u accuracy is numerical precision

plt.contourf(np.log10(np.abs(f1['jELABmat'][:,:]-f2['jELABmat'][:,:])),cmap=cm.coolwarm)


plt.colorbar()

"""
plt.figure()
f = h5py.File("gfunc/Gfunc_n1_-1_n2_2.200.h5", 'r')
for key in f.keys():
    plt.plot(f[key],color="black")
"""
