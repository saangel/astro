import pyfits as p
import matplotlib.pyplot as pl
from matplotlib.pylab import matshow
import matplotlib.mlab as mlab
from scipy.stats import norm
import numpy as n
from astroML.plotting import hist
from lmfit.models import SkewedGaussianModel
from sklearn.neighbors import KernelDensity as KDE
from astropy.stats import mad_std
import os
from astropy.utils.console import ProgressBar

def convert_catalog(cat_table):
	'''
	This function selects stars from a size-mag catalogue
	Inputs: 
	cat_table: directory of catalogue of all sources on an image. Must contain FWHM_IMAGE or FLUX_RADIUS; MAG_APER or FLUX_APER; FLAGS; VIGNET; etc
	cat_4PSFEx_table: FITS_LDAC catalogue ready to be processed by PSFEx.
	'''
	path=os.getcwd()
	direc=cat_table.split('/')[-1]
	if not os.path.exists(direc):
		os.makedirs(direc)
	os.chdir(direc)
	hdu=p.open(cat_table)
	data=hdu[2].data
	reff=data['FLUX_RADIUS']
	flux=data['FLUX_APER']
	mags=30-2.5*n.log10(flux)
	flags=data['FLAGS']
	mask = n.where((flux > 0) & (flags < 4) & (reff>0))[0]
	reff1=reff[mask]
	mags1=mags[mask]
	s=n.argsort(mags1)
	mags1=mags1[s]
	reff1=reff1[s]
	perc=0.1
	i=1
	medians=[]
	mad_stds=[]
	gammas=[]
	chisq=[]
	centers=[]
	sigmas=[]
	while perc<5:
		mask2=n.where(mags1<n.percentile(mags1,perc))[0]
		reff2=reff1[mask2]
		mags2=mags1[mask2]
		# N,bins,patches=hist(reff2,bins='scotts',label='clean, n_obj='+str(len(reff2)),normed=True,histtype='step')
		# hist(reff1,bins=bins,label='not so clean, n_obj='+str(len(reff1)),normed=True,histtype='step')
		# hist(reff,bins=bins,label='not clean at all, n_obj='+str(len(reff)),normed=True,histtype='step')
		# pl.legend()
		# pl.xlim(bins[0],bins[-1])
		# pl.savefig('hist_'+str(i)+'.png')
		# pl.clf()
		# print 'done hist'
		medians.append(n.median(reff2))
		mad_stds.append(mad_std(reff2))
		N,bins,patches=hist(reff2,bins='scotts',histtype='stepfilled',color='g',normed=False)
		pl.axvline(n.median(reff2),color='red')
		pl.axvline(n.median(reff2)+mad_std(reff2),color='black')
		pl.axvline(n.median(reff2)-mad_std(reff2),color='black')
		pl.axvline(n.median(reff2)+n.std(reff2),color='grey')
		pl.axvline(n.median(reff2)-n.std(reff2),color='grey')
		X=reff2[:,n.newaxis]
		X_plot=n.linspace(min(X)[0],max(X)[0],10000)[:,n.newaxis]
		kde=KDE(kernel='gaussian',bandwidth=min(mad_stds)).fit(X)
		ld=kde.score_samples(X_plot)
		pl.plot(X_plot[:,0],n.exp(ld)*min(mad_stds)*len(reff2),lw=3)
#		pl.xscale('log')
#		pl.yscale('log')
		'''
		model = SkewedGaussianModel()
		x = n.array([0.5 * (bins[i] + bins[i+1]) for i in xrange(len(bins)-1)])
		pars = model.guess(N, x=x)
		result=model.fit(N,pars,x=x)
		print(result.fit_report())
		pl.plot(x, result.best_fit,'k--',lw=3) 
		'''
		pl.savefig('hist_clean_log_'+str(i)+'.png')
		pl.clf()
		mask3=n.where((reff2 >= n.median(reff2)-mad_std(reff2)) & (reff2 <= n.median(reff2)+mad_std(reff2)))[0]
		reff3=reff2[mask3]
		N,bins,patches=hist(reff3,bins='scotts',histtype='stepfilled',color='g',alpha=.5,normed=True)
		x = n.array([0.5 * (bins[j] + bins[j+1]) for j in xrange(len(bins)-1)])
		#X=reff3[:,n.newaxis]
		x_plot=n.linspace(min(bins),max(bins),100)
		#kde=KDE(kernel='epanechnikov',bandwidth=min(mad_stds)).fit(X)
		#ld=kde.score_samples(X_plot)
		#pl.plot(X_plot[:,0],n.exp(ld),lw=3)
		model = SkewedGaussianModel()
		pars = model.guess(N, x=x)
		result=model.fit(N,pars,x=x)
		gammas.append(result.params['gamma'].value)
		chisq.append(result.redchi)
		centers.append(result.params['center'].value)
		sigmas.append(result.params['sigma'].value)
		#print result.fit_report()
		pl.plot(x_plot, result.eval(x=x_plot),'k--',lw=3) 
#		pl.axvline(result.params['center'].value,color='red')
#		pl.axvline(result.params['center'].value+result.params['sigma'].value,color='grey')
#		pl.axvline(result.params['center'].value-result.params['sigma'].value,color='grey')
		pl.axvline(result.params['center'].value,color='red')
		pl.axvline(result.params['center'].value+result.params['sigma'].value,color='black')
		pl.axvline(result.params['center'].value-result.params['sigma'].value,color='black')
		mu,sigma=norm.fit(reff3,loc=max(result.eval(x=x_plot)),scale=result.params['sigma'].value)
		pdf=norm.pdf(x_plot,loc=mu,scale=sigma)
		pl.plot(x_plot,pdf,color='yellow')
		pl.savefig('hist_stellarseq_'+str(i)+'.png')
		pl.clf()
		print 'done hist stellarseq'
	
		pl.plot(reff,mags,'b.',label='not clean at all',alpha=.5)
		#pl.plot(reff1,mags1,'k.',label='not so clean',alpha=.5)
		pl.plot(reff2,mags2,'r.',label='clean, mag_lims = '+str(n.min(mags2))+', '+str(n.max(mags2)),alpha=.5)
		pl.axvline(n.median(reff2),color='green')
		pl.axhline(n.max(mags2),color='black')
		pl.xlim([0,20])
		pl.ylim([n.percentile(mags1,99),n.percentile(mags1,0)])
		pl.legend()
		pl.savefig('magsize_'+str(i)+'.png')
		pl.clf()
		print 'done all'
		perc=perc+0.05
		print 'end loop '+str(i)
		i=i+1
	pl.clf()
	f, axes = pl.subplots(2,3,sharex=True)
	names=('median','mad_std','gamma','redchi','center','sigma')
	for ax,ind,name in zip(axes.flat,(medians,mad_stds,gammas,chisq,centers,sigmas),names):
		ax.plot(ind,'o')
		ax.set_title(name)
#	plt.show()
#	pl.plot(medians/n.max(medians),'o',label='median')
#	pl.plot(mad_stds/n.max(mad_stds),'o',label='mad_std')
#	pl.plot(gammas/n.max(gammas),'o',label='gammas')
#	pl.plot(chisq/n.max(chisq),'o',label='chisqr')
#	pl.legend(loc='best')
	pl.savefig('param_variations.png')
	pl.clf()
	os.chdir(path)
	'''
	vignets2=vignets1[idx]
	print len(vignets1), len(vignets2)
	samp=sorted(n.random.randint(len(idx),size=25))
	f,axes=pl.subplots(5,5,sharex=True,sharey=True)
	for ax,s in zip(axes.flat,samp):
		im=ax.imshow(vignets2[s])
	axes[0,0].set_xlim(0,36)
	axes[0,0].set_ylim(0,36)
	axes[0,0].axis('equal')
	f.subplots_adjust(right=.8)
	cbar_ax=f.add_axes([0.85,.15,.05,.7])
	f.colorbar(im,cax=cbar_ax)
	pl.show()
	samp=sorted(n.random.randint(nobj1,size=100))
	f,axes=pl.subplots(10,10,sharex=True,sharey=True)
	for i,s in zip(range(len(axes.flat)),samp):
		im=axes.flat[i].imshow(vignets1[s])
	axes[0,0].set_xlim(0,36)
	axes[0,0].set_ylim(0,36)
	axes[0,0].axis('equal')
	f.subplots_adjust(right=.8)
	cbar_ax=f.add_axes([0.85,.15,.05,.7])
	f.colorbar(im,cax=cbar_ax)
	pl.show()
	'''
dire='/Users/simon/Documents/'
direc=[dire+'CFHTLS-WIRDS/D2/R/Step 2/Proc/4PSFEx.fits',dire+'NGVS/NGVS/NGVS+0+1/Ks/Phase2/test.cat',dire+'NGVS/NGVS/NGVS+0+1/r/NGVS+0+1R.cat']
ProgressBar.map(convert_catalog,direc,multiprocess=True)
