import pyfits as p
import matplotlib.pyplot as pl
import glob as g
import subprocess as sp
import os
import numpy as np

filelist=sorted(g.glob("your*fits"))

command="#!/bin/csh \ndaospec << DONE >! logfile"
command2="\n\nDONE"
success=[]
vels=[]
verrs=[]
fwhms=[]
index=[]
i=0
for f in filelist:
	t=open('daoscript.sh','w')
	s=command+"_"+f[:-4]+"dat\n\n"+f+command2
	t.write(s)
	t.close()
	print "attempting to analyse "+f
	sp.call(['chmod','+X','daoscript.sh'])
	sp.call(['csh','daoscript.sh'])
	try:
		if os.stat(f[:-4]+"daospec").st_size > 0:
			print "File non empty!"
			success.append(f)
			fi=open(f[:-4]+"daospec")
			lines=fi.readlines()
			info=lines[0].split()
			vels.append(info[3])
			verrs.append(info[6])
			fwhms.append(info[9][:-1])
			fi.close()
			index.append(i)
		else:
			print "File empty :-("
	except OSError:
		print "File doesn't exist (?)"
	i+=1
vels=np.array(vels).astype(np.float)
verrs=np.array(verrs).astype(np.float)
fwhms=np.array(fwhms).astype(np.float)
mask=~np.isnan(vels)
