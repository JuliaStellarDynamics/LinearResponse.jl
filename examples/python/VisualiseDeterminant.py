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


model = 'IsochroneE'
dfname = 'isotropic'


model = 'PlummerE'
dfname = 'roiinf'
dfname = 'isotropic'
lharmonic = 1
n1max = 10
rb = 5.0
Ku = 202
Kv = 200
modedir = "/Users/mpetersen/Code/LinearResponse.jl/examples/{}/xifunc/".format(model)
filename = modedir+"Determinant_{}_df_{}_l_{}_n1_{}_rb_{}_Ku_{}_Kv_{}.h5".format(model,dfname,lharmonic,n1max,rb,Ku,Kv)

# need to specify these by hand: should be a part of WriteDeterminant
nx,ny = 51,50

omgarr,etaarr,detarr = read_det_file(filename,nx,ny)


fig = plt.figure(figsize=(4.5,3.5),facecolor='white')

fig = plt.gcf()
xmin = 0.17
ymin = 0.13
dx = 0.65
dy = 0.83

ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy])
ax2 = fig.add_axes([xmin+1*dx+0.02,ymin+0*dy,0.03,dy])

axlist = [ax1]


ax1.contourf(omgarr,etaarr,np.log10(np.abs(detarr)**2),np.linspace(-13,-3,100),cmap=cm.magma)


# this needs a colourbar
ax1.axis([0.,0.04,-0.01,0.05])
ax1.tick_params(axis="both",direction="in",which="both")
ax1.set_xlabel('Re[$\omega$]')
ax1.set_ylabel('Im[$\omega$]')


cmap = cm.magma
cmap = cmap; norm = mpl.colors.Normalize(vmin=-13., vmax=-3.)
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm)
cb1.set_label('log$_{10}$($\epsilon_{\ell=1}(\omega)$)')
cb1.set_ticks([-12.,-8.,-4.])
cb1.ax.minorticks_off()




plt.tight_layout()
plt.savefig('{}_n1max10.png'.format(model),dpi=300)


"""

using OrbitalElements

bc, M = 1.,1.
ψ(r::Float64)::Float64   = OrbitalElements.ψPlummer(r,bc,M,G)
dψ(r::Float64)::Float64  = OrbitalElements.dψPlummer(r,bc,M,G)
d2ψ(r::Float64)::Float64 = OrbitalElements.d2ψPlummer(r,bc,M,G)
Ω₀ = OrbitalElements.Ω₀Plummer(bc,M,G)

a,e = .1,0.1
n1,n2 = -1,2
o1,o2 = OrbitalElements.ComputeFrequenciesAE(ψ,dψ,d2ψ,a,e)
ndotOmega = n1*o1 + n2*o2
E,L = OrbitalElements.ELFromAE(ψ,dψ,a,e)
OrbitalElements.plummer_ISO_ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega::Float64)

OrbitalElements.plummer_ROI_ndFdJ(n1::Int64,n2::Int64,E::Float64,L::Float64,ndotOmega,1.,1.,1.,1.e5)


"""
