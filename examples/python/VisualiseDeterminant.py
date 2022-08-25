"""
this is an example script showing how to visualise determinant files.

@IMPROVE: make it so the determinant file records nx,ny -- rather than being user-specified.


"""

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



def read_det_file(filename,nx,ny):
    """custom reader for determinant files"""
    g = h5py.File(filename, 'r')
    omgarr = np.array(g['omega']).reshape(nx,ny)
    etaarr = np.array(g['eta']).reshape(nx,ny)
    detarr = np.array(g['det']).reshape(nx,ny)
    g.close()
    return omgarr,etaarr,detarr


model = 'IsochroneA'
modedir = "examples/{}/xifunc/".format(model)
filename = modedir+"Determinant_{}_df_isotropic_l_1_n1_10_rb_20.0.200.h5".format(model)

# need to specify these by hand: should be a part of WriteDeterminant
nx,ny = 51,50

omgarr,etaarr,detarr = read_det_file(filename,nx,ny)

plt.figure(figsize=(5,3),facecolor='white')

plt.contourf(omgarr,etaarr,np.log10(np.abs(detarr)**2),np.linspace(-13,-3,100),cmap=cm.magma)#,np.linspace(-10,0,100),cmap=cm.inferno)

plt.colorbar(label='log$_{10}$($\epsilon_{\ell=1}(\omega)$)')
plt.xlabel('Re[$\omega$]')
plt.ylabel('Im[$\omega$]')
plt.title('Isotropic Isochrone damped mode test')

plt.tight_layout()
plt.savefig('{}_n1max10.png'.format(model),dpi=300)
