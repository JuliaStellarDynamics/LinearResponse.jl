"""
InspectLinearResponse.py

a simple Python script to read output files from LinearResponse

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
plt.ion()

# HDF5 reader for Python
import h5py


fig = plt.figure(figsize=(3.5,3.0),facecolor='white')

fig = plt.gcf()
xmin = 0.2
ymin = 0.15
dx = 0.75
dy = 0.8

ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy])

axlist = [ax1]

omgfac = 2.0
nruns = 5

for imode,modepoint in enumerate(['_b_nrad100','_a_nrad100','_c_nrad100','_d_nrad25','_e_nrad25','_f_nrad100']):
    A = np.genfromtxt("/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/mode"+modepoint+".txt")
    ravals,mode1,det1 = A[:,0],A[:,2],np.abs(A[:,3])
    ax1.scatter(ravals[det1<1.],mode1[det1<1.]/omgfac,facecolor=cm.viridis(imode/nruns,1.),edgecolor='none',s=5.,zorder=-imode)

ax1.plot([0.,10.],[0.,0.],color='gray',lw=0.5,linestyle='dashed')

for ax in axlist:
    ax.axis([0.75,1.05,-0.005,0.1/omgfac])
    ax.tick_params(axis="both",direction="in",which="both")

xmax = 1.01
ax1.text(xmax,0.036,'Run A',color=cm.viridis(1./nruns,1.),ha='right',va='top')
ax1.text(xmax,0.033,'Run B',color=cm.viridis(0./nruns,1.),ha='right',va='top')
ax1.text(xmax,0.030,'Run C',color=cm.viridis(2./nruns,1.),ha='right',va='top')
ax1.text(xmax,0.027,'Run D',color=cm.viridis(3./nruns,1.),ha='right',va='top')
ax1.text(xmax,0.024,'Run E',color=cm.viridis(4./nruns,1.),ha='right',va='top')
ax1.text(xmax,0.021,'Run F',color=cm.viridis(4./nruns,1.),ha='right',va='top')
ax1.set_xlabel('anisotropy radius $r_a$')
ax1.set_ylabel('growth rate $\gamma$')

plt.savefig('PlummerEDemonstration.png',dpi=300)
"""
convert PlummerEDemonstration.png PlummerEDemonstration.pdf
pdftops -eps -level3 PlummerEDemonstration.pdf PlummerEDemonstration.eps
rm PlummerEDemonstration.pdf
"""



fig = plt.figure(figsize=(4.5,3.5),facecolor='white')

fig = plt.gcf()
xmin = 0.17
ymin = 0.13
dx = 0.65
dy = 0.83

ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy])
ax2 = fig.add_axes([xmin+1*dx+0.02,ymin+0*dy,0.03,dy])

axlist = [ax1]


minra = 0.75

startinggammas = ["1.0e-5",0.01,0.02,0.05]
setvals = [[202,200,1.0]] # a
setvals = [[202,200,5.0]] # b
#setvals = [[202,200,10.0]] # c
#setvals = [[51,50,5.0]] # e
#setvals = [[202,200,20.0]] # f

cmin,cmax = 0.75, 1.05
for iset,setval in enumerate(setvals):
    Ku,Kv,rb = setval
    for gammastart in startinggammas:
        for ra in range(1,50,1):#120):
            raval = np.round((ra)*0.005 + minra,3)
            modefile = "/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/xifunc/ModeShape_{4}_PlummerE_df_roi{0}_l_2_n1_1_rb_{1}_Ku_{2}_Kv_{3}.h5".format(raval,rb,Ku,Kv,gammastart)
            #print(modefile)
            try:
                f = h5py.File(modefile, 'r')
                print(modefile)
                dvals = np.real(f['ModePotentialShape'][:])
                peakdens = np.nanargmax(np.abs(dvals))
                if np.nanmin(dvals/dvals[peakdens])<-0.1:
                    if raval<0.85:
                        _ = ax1.plot(f['ModeRadius'][:],0.5*dvals/dvals[peakdens],color=cm.viridis(((raval-cmin)/(cmax-cmin))),linestyle='dashed')
                else:
                    _ = ax1.plot(f['ModeRadius'][:],dvals/dvals[peakdens],color=cm.viridis(((raval-cmin)/(cmax-cmin))))
            except:
                pass


# this needs a colourbar
ax1.axis([0.,12.0,-0.35,1.05])
ax1.tick_params(axis="both",direction="in",which="both")
ax1.set_xlabel('cluster radius $r$')
ax1.set_ylabel('potential fluctuation (normalised)')


cmap = cm.viridis
cmap = cmap; norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm)
cb1.set_label('$r_a$')
#cb1.set_ticks([-1.,0.,1.])
cb1.ax.minorticks_off()


plt.savefig('PlummerECurveDemonstration.png',dpi=300)
"""
convert PlummerECurveDemonstration.png PlummerEDemonstration.pdf
pdftops -eps -level3 PlummerEDemonstration.pdf PlummerECurveDemonstration.eps
rm PlummerEDemonstration.pdf
"""





fig = plt.figure(figsize=(8.5,2.6),facecolor='white')

fig = plt.gcf()
xmin = 0.08
ymin = 0.165
dx = 0.34
dy = 0.78
xbuf = 0.13

ax1 = fig.add_axes([xmin+0*dx          ,ymin+0*dy,dx,dy])
ax2 = fig.add_axes([xmin+1*dx+xbuf     ,ymin+0*dy,dx,dy])
ax3 = fig.add_axes([xmin+2*dx+xbuf+0.02,ymin+0*dy,0.01,dy])


axlist = [ax1,ax2]

omgfac = 1.0
nruns = 5

modeshape = ['x','^','o','+','v','s']
modelabel = ['A','B','C','D','E','F']

cmin,cmax = 0.75, 1.05

for imode,modepoint in enumerate(['_a_nrad100','_b_nrad100','_c_nrad100','_d_nrad25','_e_nrad25','_f_nrad100']):
    A = np.genfromtxt("/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/mode"+modepoint+".txt")
    ravals,mode1,det1 = A[:,0],A[:,2],np.abs(A[:,3])
    if imode==1:
        ax1.scatter(ravals[det1<1.],mode1[det1<1.],marker=modeshape[imode],facecolor=cm.viridis((ravals[det1<1.]-cmin)/(cmax-cmin),1.),edgecolor='none',s=15.,zorder=-imode)
    else:
        ax1.scatter(ravals[det1<1.],mode1[det1<1.],marker=modeshape[imode],facecolor='black',edgecolor='none',s=5.,zorder=10)


ax1.plot([0.,10.],[0.,0.],color='gray',lw=0.5,linestyle='dashed')


xmax = 0.98
ymax = 0.96
dy = 0.06
for imode in range(0,6):
    if imode==1:
        ax1.text(xmax-0.02,ymax-imode*dy,'Run {}'.format(modelabel[imode]),color='black',ha='right',va='center',transform=ax1.transAxes)
        ax1.scatter([xmax],[ymax-imode*dy],marker=modeshape[imode],color='black',transform=ax1.transAxes,s=10.)
    else:
        ax1.text(xmax-0.02,ymax-imode*dy,'Run {}'.format(modelabel[imode]),color='grey',ha='right',va='center',transform=ax1.transAxes)
        ax1.scatter([xmax],[ymax-imode*dy],marker=modeshape[imode],color='grey',transform=ax1.transAxes,s=10.)


ax1.axis([0.75,1.05,-0.005,0.12])
ax1.set_xlabel('anisotropy radius $r_a$')
ax1.set_ylabel('growth rate $\gamma$')


minra = 0.75
startinggammas = ["1.0e-5",0.01,0.02,0.05]
setvals = [[202,200,5.0]] # b

for iset,setval in enumerate(setvals):
    Ku,Kv,rb = setval
    for gammastart in startinggammas:
        for ra in range(1,50,1):#120):
            raval = np.round((ra)*0.005 + minra,3)
            modefile = "/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/xifunc/ModeShape_{4}_PlummerE_df_roi{0}_l_2_n1_1_rb_{1}_Ku_{2}_Kv_{3}.h5".format(raval,rb,Ku,Kv,gammastart)
            #print(modefile)
            try:
                f = h5py.File(modefile, 'r')
                print(modefile)
                dvals = np.real(f['ModePotentialShape'][:])
                peakdens = np.nanargmax(np.abs(dvals))
                if np.nanmin(dvals/dvals[peakdens])<-0.1:
                    if raval<0.85:
                        _ = ax2.plot(f['ModeRadius'][:],0.5*dvals/dvals[peakdens],color=cm.viridis(((raval-cmin)/(cmax-cmin))),linestyle='dashed')
                else:
                    _ = ax2.plot(f['ModeRadius'][:],dvals/dvals[peakdens],color=cm.viridis(((raval-cmin)/(cmax-cmin))))
            except:
                pass


# this needs a colourbar
ax2.axis([0.,12.0,-0.35,1.05])
ax2.set_xlabel('model radius $r$')
ax2.set_ylabel('potential fluctuation (normalised)')

for ax in axlist:
  ax.tick_params(axis="both",direction="in",which="both")

cmap = cm.viridis
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,norm=norm)
cb1.set_label('$r_a$')
#cb1.set_ticks([0.75,0.85,0.95])
cb1.ax.minorticks_off()


ax1.text(0.02,0.98,'(a)',color='black',ha='left',va='top',transform=ax1.transAxes)
ax2.text(0.98,0.98,'(b)',color='black',ha='right',va='top',transform=ax2.transAxes)

plt.savefig('PlummerECurveDemonstration.png',dpi=300)
"""
convert PlummerECurveDemonstration.png PlummerEDemonstration.pdf
pdftops -eps -level3 PlummerEDemonstration.pdf PlummerECurveDemonstration.eps
rm PlummerEDemonstration.pdf
"""



modefile = "/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/xifunc/Determinant_PlummerE_df_isotropic_l_1_n1_10_rb_1.0_Ku_202_Kv_200.h5"

modefile = "/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/xifuncDeterminant_PlummerE_df_isotropic_l_1_n1_10_rb_20.0_Ku_202_Kv_200.h5"
f = h5py.File(modefile, 'r')

plt.contourf(np.array(f['omega']).reshape(51,50),np.array(f['eta']).reshape(51,50),np.log10(np.abs(np.array(f['det'])).reshape(51,50)),np.linspace(-12.,-4.,10))

plt.xlabel('omega')
plt.xlabel('omega')



filename = '/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/gfunc/Gfunc_PlummerE_df_roi0.75_l_1_n1_-1_n2_-1_rb_4.0_Ku_202_Kv_200.h5'
filename = '/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/gfunc/Gfunc_PlummerE_df_roi0.75_l_1_n1_2_n2_-1_rb_4.0_Ku_202_Kv_200.h5'

filename = '/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/gfunc/Gfunc_PlummerE_df_roi0.75_l_1_n1_-1_n2_2_rb_4.0_Ku_202_Kv_200.h5'

filename = '/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/gfunc/Gfunc_PlummerE_df_roi0.75_l_2_n1_-1_n2_2_rb_4.0_Ku_202_Kv_200.h5'

#filename = '/Users/mpetersen/Code/LinearResponse.jl/examples/PlummerE/gfunc/Gfunc_PlummerE_df_roi0.75_l_2_n1_1_n2_-2_rb_4.0_Ku_202_Kv_200.h5'

f = h5py.File(filename, 'r')
print(f.keys())

n,w = np.polynomial.legendre.leggauss(202)
plt.plot(n,f['Gmat'][:,1,1])

z = np.polyfit(n,f['Gmat'][:,1,1], 30)
fz = np.poly1d(z)
plt.plot(n,vz(n))

z1 = np.polyder(z,1)
fz1 = np.poly1d(z1)
plt.plot(n,fz1(n))



deriv =


vmin,vmax = OrbitalElements.FindVminVmax(uval,n1,n2,dψ,d2ψ,ωmin,ωmax,Orbitalparams)
