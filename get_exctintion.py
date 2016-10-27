import astropy.coordinates as coord
import astropy.units as u
import glob as g
import numpy as np
from extinction import fitzpatrick99 as f99
import matplotlib.pyplot as plt
from astroquery.irsa_dust import IrsaDust


coords=["02h25m59.00s    -04d29m40.00s","10h00m28.00s    02d12m30.00s","14h19m27.00s    52d40m56.00s","22h15m31.00s    -17d43m56.00s"]
filter_names="ugrizJHK"
filters=[]
for f in filter_names:
    filters.append((g.glob("synth/*"+f+"*Mega*")+g.glob("synth/"+f+"wircam*"))[0])
wl=np.empty(0)
filts=[]
wls=[]
ax1=plt.subplot2grid((2,1),(0,0))
ax2=plt.subplot2grid((2,1),(1,0),sharex=ax1)
spec=0.5
for f in filters:
    wlr,eff=np.loadtxt(f,unpack=True)
    wls.append(wlr)
    wl=np.concatenate((wl,wlr))
    filts.append(eff)
    ax1.fill(wlr,eff,label=f.split('.')[0],alpha=.5,edgecolor="none")
    ax1.axhline(spec,color="black",lw=3,alpha=.5)
#    ax1.set_xlabel(r"$\lambda$ in $\AA$")
    ax1.set_ylabel("Throughput")
    ax1.axes.get_xaxis().set_visible(False)
wl=np.sort(wl)

corrections=np.empty((len(filters),len(coords)))
mags_notred=np.empty(len(filters))
mags_red=np.empty((len(filters),len(coords)))
alambdas=[ [[] for _ in coords] for _ in filts]

for i,c in enumerate(coords):
  C = coord.SkyCoord(c,frame="fk5")
  table=IrsaDust.get_query_table(c, radius=2.0 * u.deg)
  eb_v=table["ext SandF mean"]
  #print eb_v.data[0]
  al_plot=f99(wl,eb_v.data[0]*3.1)
  for j,f in enumerate(filts):
      alambdas[j][i]=f99(wls[j],eb_v.data[0]*3.1)
  ax2.plot(wl,al_plot,label="D"+str(i+1))
  ax2.set_xlabel(r"$\lambda$ in $\rm \AA$")
  ax2.set_ylabel("Extinction in magnitudes")
alambdas=np.array(alambdas)

for j,f in enumerate(filts):
    diffs=np.gradient(wls[j])
    flux=sum(wls[j]*spec*f*diffs) #integration
    norm=sum(f*diffs/wls[j]) #normalisation following GALEXEV docsself.
    for k,c in enumerate(coords):
        tau=alambdas[j,k]/2.5/np.log10(np.e)
        spec_red=spec*np.exp(-tau)
        flux_red=sum(wls[j]*spec_red*f*diffs)
        mags_red[j,k]=-2.5*np.log10(flux_red/norm)
    mags_notred[j]=-2.5*np.log10(flux/norm)

print filter_names
derred=[]
for j in range(len(coords)):
    print mags_notred-mags_red[:,j]
    derred.append(mags_notred-mags_red[:,j])
np.savetxt("derredening_indexes.txt",derred)
plt.legend(loc="best")
title=" ".join([filters[i].split(".")[0][6:]+" = "+str(m)[:8] for i,m in enumerate(mags_notred)])
title=title+"\n"
title=title+" ".join([filters[i].split(".")[0][6:]+" = "+str(m)[:8] for i,m in enumerate(mags_red[:,1])])
#plt.title(title)
plt.show()
