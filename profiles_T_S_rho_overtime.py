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

pii=3.14159

fmesh=netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
z = np.zeros(depths.shape)
z[0] = -0.5 * depths[0]
z[1:] = -0.5 * (depths[0:-1] + depths[1:])
#nz = len(z)

#idx=67250-1
idx = np.argmin( (latCell-(-1.323514))**2 + (lonCell-5.672896)**2)  #122901-1
#idx=198673-1
#idx=210384-1 # filchner sill
print "idx=",idx

maxLevelCell=fmesh.variables['maxLevelCell'][idx]
nz = maxLevelCell



# plot location
fig = plt.figure(2, facecolor='w')
idx2 = np.nonzero(latCell<-60.0/180.0*pii)[0]
plt.plot(yCell[idx2], xCell[idx2], 'k.')
plt.plot(yCell[idx], xCell[idx], 'r.')

#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'


years = np.arange(50,113,1)
months = np.arange(1,13,1)
nt = len(years)*len(months)
times = np.zeros((nt,))


Tdata = np.zeros((nz, nt))
Sdata = np.zeros((nz, nt))
rhodata = np.zeros((nz, nt))
vdata = np.zeros((nz, nt))

t=0
for yr in years:
 for mo in months:
   print "yr=",yr, "mo=", mo
   times[t] = yr+(mo-1.0)/12.0

   f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   T=f.variables['timeMonthly_avg_activeTracers_temperature']
   S=f.variables['timeMonthly_avg_activeTracers_salinity']
   rho = f.variables['timeMonthly_avg_potentialDensity']
   v  = f.variables['timeMonthly_avg_velocityMeridional']
   
   # get data
   Tdata[:,t] = T[0,idx,:nz]
   Sdata[:,t] = S[0,idx,:nz]
   rhodata[:,t] = rho[0,idx,:nz]
   vdata[:,t] = v[0,idx,:nz]

   f.close()
   t += 1
  
fig = plt.figure(1, facecolor='w')
nrow=4
ncol=1

maxDepth=-200.0

axT = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('year')
plt.ylabel('temperature (deg. C)')
plt.pcolor(times, z[:nz], Tdata, vmin=-2.1, vmax=1.6, cmap='nipy_spectral')#, vmin=-2.0, vmax=0.0)
#axT.set_ylim((maxDepth, 0))
plt.colorbar()
#plt.contour(times, z[:nz], Tdata, [-1.8])

axS = fig.add_subplot(nrow, ncol, 2)
plt.xlabel('year')
plt.ylabel('salinity (psu)')
plt.pcolor(times, z[:nz], Sdata, vmin=33.5,vmax=34.6, cmap='nipy_spectral')#, vmin=31.0, vmax=34.7)
#axS.set_ylim((maxDepth, 0))
plt.colorbar()

axrho = fig.add_subplot(nrow, ncol, 3)
plt.xlabel('year')
plt.ylabel('potential density (kg/m3)')
plt.pcolor(times, z[:nz], rhodata, vmin=1027.2,vmax=1027.7, cmap='nipy_spectral')#, vmin=31.0, vmax=34.7)
#axS.set_ylim((maxDepth, 0))
plt.colorbar()

axv = fig.add_subplot(nrow, ncol, 4)
plt.xlabel('year')
plt.ylabel('northward velo (m/s)')
plt.pcolor(times, z[:nz], vdata, cmap='seismic', vmin=-0.05, vmax=0.05)
#axS.set_ylim((maxDepth, 0))
plt.colorbar()



plt.show()
