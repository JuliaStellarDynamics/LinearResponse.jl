
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
f = h5py.File('wmat/wmat_l_2_n1_-1_n2_-2.h5', 'r')


f = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/test/Wmat_Plummer_l_2_n1_0_n2_2_rb_3.0_Ku_12_Kv_20_Kw_20.h5', 'r')

plt.contourf(f['wmat'][:,:,3],cmap=cm.coolwarm)
#plt.contourf(np.log10(f['amat'][:,:]),64,cmap=cm.coolwarm)
#plt.contourf((f['emat'][:,:]),64,cmap=cm.coolwarm)

plt.colorbar()

f = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/test/Wmat_Plummer_l_2_n1_-1_n2_2_rb_3.0_Ku_52_Kv_20_Kw_20.h5', 'r')

alpha = f['Omgmat'][:,:,0].flatten()
beta = (f['Omgmat'][:,:,1]/f['Omgmat'][:,:,0]).flatten()
cmin,cmax = np.nanmin(f['wmat'][:,:,3]),np.nanmax(f['wmat'][:,:,3])
plt.scatter(alpha,beta,c=cm.coolwarm((f['wmat'][:,:,3].flatten()-cmin)/(cmax-cmin)))


f = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/test/Wmat_Plummer_l_2_n1_0_n2_2_rb_3.0_Ku_52_Kv_20_Kw_20.h5', 'r')
alpha = f['Omgmat'][:,:,0].flatten()
beta = (f['Omgmat'][:,:,1]/f['Omgmat'][:,:,0]).flatten()
plt.scatter(alpha,beta,c=cm.coolwarm((f['UVmat'][:,:,0].flatten()+1)/(2)))


f = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/test/Wmat_Plummer_l_2_n1_0_n2_2_rb_3.0_Ku_52_Kv_20_Kw_20.h5', 'r')
alpha = f['Omgmat'][:,:,0].flatten()
beta = (f['Omgmat'][:,:,1]/f['Omgmat'][:,:,0]).flatten()
plt.scatter(alpha,beta,c=cm.coolwarm((f['UVmat'][:,:,0].flatten()+1)/(2)))

f = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/test/Wmat_Plummer_l_2_n1_4_n2_2_rb_3.0_Ku_52_Kv_20_Kw_20.h5', 'r')
alpha = f['Omgmat'][:,:,0].flatten()
beta = (f['Omgmat'][:,:,1]/f['Omgmat'][:,:,0]).flatten()
plt.scatter(alpha,beta,c=cm.coolwarm((f['UVmat'][:,:,0].flatten()+1)/(2)))


f = h5py.File('/Users/mpetersen/Code/LinearResponse.jl/test/Wmat_Plummer_l_2_n1_-3_n2_2_rb_3.0_Ku_52_Kv_20_Kw_20.h5', 'r')
alpha = f['Omgmat'][:,:,0].flatten()
beta = (f['Omgmat'][:,:,1]/f['Omgmat'][:,:,0]).flatten()
plt.scatter(alpha,beta,c=cm.coolwarm((f['UVmat'][:,:,0].flatten()+1)/(2)))






plt.figure()
f = h5py.File("/Users/mpetersen/Code/LinearResponse.jl/test/Gfunc_Plummer_df_roi0.75_l_2_n1_-1_n2_2_rb_3.0_Ku_52_Kv_20.h5", 'r')
plt.plot(uvals,f['Gmat'][:,0,0]/((uvals+1)*(uvals-1)),color="black")


f = h5py.File("/Users/mpetersen/Code/LinearResponse.jl/test/Gfunc_Plummer_df_roi0.75_l_2_n1_-3_n2_2_rb_3.0_Ku_52_Kv_20.h5", 'r')
plt.plot(f['Gmat'][:,0,0],color="black")
