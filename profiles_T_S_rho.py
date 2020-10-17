#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
'''

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm
from math import pi

#run = 'ISMF'
run = 'ISMF-noEAIS'
startyr = 82
endyr = 82

savepath='/global/homes/c/cbegeman/weddell_output/'

fmesh=netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
runname = ['ISMF',
           'ISMF-noEAIS']
runpath = ['/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
           '/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run']

latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
z = np.zeros(depths.shape)
z[0] = -0.5 * depths[0]
z[1:] = -0.5 * (depths[0:-1] + depths[1:])

deg2rad  = pi/180.
#latplt = -1.323514
#lonplt = 5.672896
latS = 73.
lonW = 32.
locname = str(latS) + 'S' + str(lonW) + 'W'
latplt = -1.*latS*deg2rad
lonplt = (360.-lonW)*deg2rad
#idx=67250-1
idx = np.argmin( (latCell-latplt)**2 + (lonCell-lonplt)**2)  #122901-1

maxDepth = -500.

years = np.arange(startyr,endyr+1,1)
months = np.arange(1,13,1)
nt = len(years)*len(months)
times = np.zeros((nt,))
colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, 13)]
#colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, nt)]

fig = plt.figure(1, facecolor='w')
nrow=1
ncol=3

axT = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('temperature (deg. C)')
plt.ylabel('depth (m)')
plt.grid()

axS = fig.add_subplot(nrow, ncol, 2, sharey=axT)
plt.xlabel('salinity (psu)')
#plt.ylabel('depth (m)')
axS.set_yticklabels([])
plt.grid()

axrho = fig.add_subplot(nrow, ncol, 3, sharey=axT)
plt.xlabel('potential density (kg/m^3)')
#plt.ylabel('depth (m)')
axrho.set_yticklabels([])
plt.grid()

lineStyle = ['_',':','-.','--']
t = 0
for i,yr in enumerate(years):
 for j,mo in enumerate(months):
   #for i in range(len(years)):
   #yr = years[i]
   #c = colors[t]
   c = colors[j]
   
   datestr = '{0:04d}-{1:02d}'.format(yr, mo)
   filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
   f = netCDF4.Dataset(filename, 'r')
   #f=netCDF4.Dataset(path+'/'+'mpaso.hist.am.timeSeriesStatsMonthly.{0:04d}-01-01.nc'.format(yr),'r')
   T=f.variables['timeMonthly_avg_activeTracers_temperature']
   S=f.variables['timeMonthly_avg_activeTracers_salinity']
   rho = f.variables['timeMonthly_avg_potentialDensity']

   axT.plot(T[0,idx,:], z, label="yr{0:04d}".format(yr), color=c,linestyle=lineStyle[i])
   #axT.plot(-0.0575*S[0,idx,:]+0.0901-7.61e-4*z, z, 'k.')
   axS.plot(S[0,idx,:], z, label="yr{0:04d}".format(yr), color=c,linestyle=lineStyle[i])
   axrho.plot(rho[0,idx,:], z, label="yr{0:04d}".format(yr), color=c,linestyle=lineStyle[i])

   t = t+1

axT.set_xlim([-2.2,2.5])
axT.set_ylim([maxDepth,0])
axS.set_xlim([32,35])
axS.set_ylim([maxDepth,0])
axrho.set_xlim([1025.8,1028])
axrho.set_ylim([maxDepth,0])

#axT.legend()
#axS.legend()

#axT.plot([-1.8, -1.8], [-4000.0, 0.0], 'k:')

plt.savefig(savepath + run + '_profiles_' + locname + '_' + str(startyr) + '-' + str(endyr) + '.png')
plt.close()
