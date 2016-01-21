import pyfits as p
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import glob as g
from joblib import Parallel,delayed
import time as t
import sys 
import copy as cp
from scipy import interpolate

speclist=g.glob('*/*.fits') #grab all spectra
filters=g.glob('*/*.txt') # grab filters
# color=cm.jet(np.linspace(0,1,len(filters))) #assign a color to each filter, this was used when spectra was plotted
wlr0=np.arange(3000,10000,0.1) #wavelenght range for PHOENIX medium resolution spectra
#The following three lines are the correction from vacuum to air wavelenghts, following Husser+13
sigma=1e8/np.power(wlr0,2)
f=1+0.05792105/(238.0185-sigma)+0.00167917/(57.362-sigma)
wlr=wlr0/f
fluxes=[0]*len(filters)
norms=[0]*len(filters)
names=[f[0] for f in filters]

def magic(spec):
	sp=p.open(spec) #load spectra
	data=sp[0].data
	for j,f in enumerate(filters):
		wl,eff=np.loadtxt(f,unpack=True) #load filter
		ff=interpolate.interp1d(wl,eff)
		mask2=np.where((wlr<=max(wl)) & (wlr>=min(wl)))
		wlr2=wlr[mask2]
		eff2=ff(wlr2)
		wls=0.1 #PHOENIX resolution. Because of the conversion above, it's not exactly 0.1, it's more like 0.099973 something
		fluxes[j]=sum(wlr2*data[mask2]*eff2*wls) #integration
		norms[j]=sum(wlr2*eff2*wls) #normalisation following GALEXEV docs.
	s=np.argsort(names)
	filters_=np.array(filters)
	fluxes_=np.array(fluxes)
	filters_=filters_[s]
	fluxes_=fluxes_[s]
	# The following commented lines plotted the spectra
	# plt.title("u-g={:.3},g-r={:.3},r-i={:.3},i-z={:.3}".format(-2.5*np.log10(fluxes_[3]/fluxes_[0]),\
	# 	-2.5*np.log10(fluxes_[0]/fluxes_[2]),-2.5*np.log10(fluxes_[2]/fluxes_[1]),\
	# 	-2.5*np.log10(fluxes_[1]/fluxes_[4])))
	# plt.plot(wlr,data,'r-',label=spec.split('/')[1])
	# plt.xlabel('angstrom')
	# plt.ylabel('flux')
	# plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.01),fancybox=True, shadow=True)
	# plt.savefig("images/"+spec.split('/')[1]+".png",dpi=300)
	# plt.clf()
	i=0
	# Some of the PHOENIX spectra also have different values for  alpha elements abundance, this changes the name of the files
	if 'Alpha' in spec:
		i=12
	temp=int(spec[i+44:i+49])
	logg=float(spec[i+50:i+54])
	feh=float(spec[36:40])
	return [-2.5*np.log10(fluxes_[3]*norms[0]/(fluxes_[0]*norms[3])),
		-2.5*np.log10(fluxes_[0]*norms[2]/(fluxes_[2]*norms[0])),
		-2.5*np.log10(fluxes_[2]*norms[1]/(fluxes_[1]*norms[2])),
	 	-2.5*np.log10(fluxes_[1]*norms[4]/(fluxes_[4]*norms[1])),temp,logg,feh]
skips=1 # line added to do test quickly in a subset of the files
a=Parallel(n_jobs=-1,verbose=20)(delayed(magic)(s) for s in speclist[::skips]) # parallel running

# load the actual CFHTLS data
h=p.open('../Data/megacat.cat')
coords=['ALPHA_J2000_1','DELTA_J2000_1']
filters=['u','g','r','i','z']
mags=np.zeros((len(h)-1,len(filters)),dtype=object)
pos=np.zeros(len(h)-1,dtype=object)
mask=np.zeros(len(filters),dtype=object)

# i is field, j is filter
for i in range(len(mags)):
	for j in range(len(filters)):
		mags[i,j]=h[i+1].data[filters[j]]
		pos[i]=zip(h[i+1].data[coords[0]],h[i+1].data[coords[1]])
		mask[i]=(mags[i,0]<30) & (mags[i,1]<30) & (mags[i,2]<30) & (mags[i,3]<30) & (mags[i,4]<30)
pos=np.array(pos)
for i in range(len(mags)):
	pos[i]=np.array(pos[i])
for i in range(len(mags)):
	for j in range(len(filters)):
		mags[i,j]=mags[i,j][mask[i]]
	pos[i]=pos[i][mask[i]]

# reddening for each field and filter
Av=np.array([[0.113,0.088,0.061,0.045,0.034],[0.079,0.061,0.042,0.032,0.023],[0.035,0.038,0.019,0.014,0.011],[0.106,0.083,0.057,0.042,0.032]])
for i,f in enumerate(mags):
	for j,m in enumerate(f):
		mags[i,j]=m-Av[i,j]

a=np.array(a) # so I can access quickly the columns
colors=['u-g','g-r','r-i','i-z']
props=['temp','logg','feh']
for i in range(3):
	plt.plot(a[:,i],a[:,i+1],'.')
	plt.xlabel(colors[i])
	plt.ylabel(colors[i+1])
	plt.savefig('cc'+str(i+1)+'.png',dpi=300)
	plt.clf()
	plt.close('all')
	for j,p in enumerate(props):
		fig, ax = plt.subplots()
		cax = ax.scatter(a[:,i],a[:,i+1],c=a[:,4+j],s=50,edgecolor='none',zorder=1)
		for k,f in enumerate(mags):
			ax.plot(f[i]-f[i+1],f[i+1]-f[i+2],marker='1',linestyle='none',label='D'+str(1+k),alpha=0.1,color='0.25',zorder=2)
		cbar = fig.colorbar(cax)
		plt.xlabel(colors[i])
		plt.ylabel(colors[i+1])
		cbar.set_label(p)
		#plt.legend(loc='best')
		plt.savefig('cc'+str(i+1)+'+'+p+'.png',dpi=300)
		plt.clf()
		plt.close('all')

plt.close('all')
# urz plot
# for j,p in enumerate(props):
# 	fig,ax=plt.subplots()
# 	cax=ax.scatter(a[:,0]+a[:,1],a[:,2]+a[:,3],c=a[:,4+j],s=50,edgecolor='none')
# 	cbar=fig.colorbar(cax)
# 	plt.savefig('urz+props+'+p+'.png')
# 	plt.clf()
# 	plt.close('all')
mets=np.unique(a[:,-1])
f, axes = plt.subplots(3,3, sharex = True,sharey=True)
for i,ax in enumerate(axes.flat):
	m=np.where((a[:,-1]==mets[i]))
	cax=ax.scatter(a[:,1][m],a[:,2][m],c=a[:,4][m],s=5*a[:,5][m],edgecolor='none',label='[Fe/H]='+str(mets[i]))
	ax.legend(loc='best')
cbar = fig.colorbar(cax)
plt.show()