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

fmesh=netCDF4.Dataset('oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
z = np.zeros(depths.shape)
z[0] = -0.5 * depths[0]
z[1:] = -0.5 * (depths[0:-1] + depths[1:])
pii=3.14159



yr=94
idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-78.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>300.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell

print "Found {} cells".format(len(idx))

fig = plt.figure(2, facecolor='w')
idx2 = np.nonzero(latCell<-60.0/180.0*pii)[0]
plt.plot(yCell[idx2], xCell[idx2], 'k.')
plt.plot(yCell[idx], xCell[idx], 'r.')

path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'
print "running year:", yr

fig = plt.figure(1, facecolor='w')
nrow=1
ncol=1

for yr in np.arange(95,102,1):
 for mo in np.arange(1,13,1):
   print "yr=",yr, "mo=", mo
   axTS = fig.add_subplot(nrow, ncol, 1)
   plt.ylabel('temperature (deg. C)')
   plt.xlabel('salinity (psu)')
   plt.grid()


   #f=netCDF4.Dataset('/usr/projects/climate/mhoffman/e3sm/cryo-campaign-data/mpaso.hist.am.timeSeriesStatsMonthly.0070-01-01.nc')
   f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   T=f.variables['timeMonthly_avg_activeTracers_temperature']
   S=f.variables['timeMonthly_avg_activeTracers_salinity']
   
   # get data
   Ts = T[0,idx,:]
   Ss = S[0,idx,:]
   sc=axTS.scatter(Ss[:], Ts[:], s=1, c=np.tile(z, (len(idx),1)), vmin=-700, vmax=0.0, cmap='twilight')
   
   axTS.set_ylim([-2.1,2.0])
   axTS.set_xlim([33.0,34.8])
   plt.colorbar(sc)
   
   plt.plot([34.0, 34.0], [-1.85, 0.2], 'k:')
   plt.plot([34.5, 34.5], [-1.86, -1.5], 'k:')
   plt.plot([34.0, 35.0], [-1.5, -1.5], 'k:')
   plt.plot([34.0, 35.0], [0.0, 0.0], 'k:')
   plt.plot([33.5, 35.0], [-0.0575*33.5+0.0901, -0.0575*35.0+0.0901], 'k:')
   plt.title("year {}, mo {}".format(yr, mo))
   plt.draw()
   plt.savefig('TS_yr{:04d}-{:02d}.png'.format(yr, mo))
   plt.clf()
   
   
   
   # velo plots
   v=f.variables['timeMonthly_avg_velocityMeridional'][0,idx,:]
   u=f.variables['timeMonthly_avg_velocityZonal'][0,idx,:]
   swsv=f.variables['timeMonthly_avg_windStressMeridional'][0,idx]
   swsu=f.variables['timeMonthly_avg_windStressZonal'][0,idx]
   
   fig = plt.figure(3, facecolor='w')
   
   axUsrf = fig.add_subplot(1, 3, 1)
   axUsrf.axis('equal')
   plt.quiver(yCell[idx], xCell[idx], u[:,0], v[:,0])
   print "d={}, z={}".format(0, z[0])
   plt.title("year={}, depth={}m".format(yr, z[0]))
   plt.title
   
   
   d=np.argmin(np.absolute(z--250.0))
   print "d={}, z={}".format(d, z[d])
   
   axU2 = fig.add_subplot(1, 3, 2)
   axU2.axis('equal')
   plt.quiver(yCell[idx], xCell[idx], u[:,d], v[:,d])
   plt.title("year={}, depth={}m".format(yr, z[d]))
   
   axWS = fig.add_subplot(1, 3, 3)
   axWS.axis('equal')
   plt.quiver(yCell[idx], xCell[idx], swsu[:], swsv[:])
   plt.title("year={}, surface stress".format(yr))
   plt.draw()
   plt.savefig('velo_yr{:04d}.png'.format(yr))
   plt.clf()

#plt.show()
