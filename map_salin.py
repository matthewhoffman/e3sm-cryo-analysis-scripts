#!/usr/bin/env python
'''
'''

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm

fmesh=netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
z = np.zeros(depths.shape)
z[0] = -0.5 * depths[0]
z[1:] = -0.5 * (depths[0:-1] + depths[1:])
pii=3.14159


size=28
fsz=(7,5)
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-78.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>300.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>310.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #eastern weddell
idx = np.nonzero( (latCell<-65.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<20.0/360.0*2*pii))[0]; size=6;  #entire weddell+
idx = np.nonzero( (latCell<-60.0/180.0*pii) )[0]; size=1; sz=(14,10); #entire SO
idx = np.nonzero( (latCell<-60.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<80.0/360.0*2*pii))[0]; size=1.4; fsz=(15,7)  # weddell to Amery 

idxSill=210384-1 # filchner sill

path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'

diffpath = False
diffpath='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'

years = np.arange(85,110,1)
months = np.arange(1,13,1)
#months = np.arange(1,2,1)
nt = len(years)*len(months)


fig = plt.figure(1, facecolor='w', figsize=fsz)
nrow=2
ncol=2

diffRange = 0.2

k100 = np.argmin(np.absolute(z--100.0))
k50 = np.argmin(np.absolute(z--50.0))
k200 = np.argmin(np.absolute(z--200.0))
k300 = np.argmin(np.absolute(z--300.0))

for yr in years:
 for mo in months:
   print "yr=",yr, "mo=", mo

   f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   S50 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k50]
   S100 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k100]
   S200 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k200]
   S300 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k300]
   f.close()

   if diffpath:
      f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(diffpath, yr, mo), 'r')
      Sd50 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k50]
      Sd100 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k100]
      Sd200 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k200]
      Sd300 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k300]
      f.close()


   ax50 = fig.add_subplot(nrow, ncol, 1)
   if diffpath:
      sc1 = ax50.scatter(yCell[idx], xCell[idx], s=size, c=(S50-Sd50), vmin=-1.0*diffRange, vmax=diffRange, cmap='RdBu_r')
   else:
      sc1 = ax50.scatter(yCell[idx], xCell[idx], s=size, c=(S50), vmin=33.5, vmax=34.6)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax50)
   ax50.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax50.axis('equal')
   ax50.set_title("{}:\n year {}, mo {}; Sal. 50m".format(path, yr, mo), fontsize=7)

   ax100 = fig.add_subplot(nrow, ncol, 2)
   if diffpath:
      sc1 = ax100.scatter(yCell[idx], xCell[idx], s=size, c=(S100-Sd100), vmin=-1.0*diffRange, vmax=diffRange, cmap='RdBu_r')
   else:
      sc1 = ax100.scatter(yCell[idx], xCell[idx], s=size, c=(S100), vmin=33.5, vmax=34.6)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax100)
   ax100.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax100.axis('equal')
   ax100.set_title("{}:\n year {}, mo {}; Sal. 100m".format(path, yr, mo), fontsize=7)

   ax200 = fig.add_subplot(nrow, ncol, 3)
   if diffpath:
      sc1 = ax200.scatter(yCell[idx], xCell[idx], s=size, c=(S200-Sd200), vmin=-1.0*diffRange, vmax=diffRange, cmap='RdBu_r')
   else:
      sc1 = ax200.scatter(yCell[idx], xCell[idx], s=size, c=(S200), vmin=33.5, vmax=34.6)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax200)
   ax200.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax200.axis('equal')
   ax200.set_title("{}:\n year {}, mo {}; Sal. 200m".format(path, yr, mo), fontsize=7)

   ax300 = fig.add_subplot(nrow, ncol, 4)
   if diffpath:
      sc1 = ax300.scatter(yCell[idx], xCell[idx], s=size, c=(S300-Sd300), vmin=-1.0*diffRange, vmax=diffRange, cmap='RdBu_r')
   else:
      sc1 = ax300.scatter(yCell[idx], xCell[idx], s=size, c=(S300), vmin=33.5, vmax=34.6)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax300)
   ax300.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax300.axis('equal')
   ax300.set_title("{}:\n year {}, mo {}; Sal. 300m".format(path, yr, mo), fontsize=7)



   plt.draw()
   plt.savefig('salin_map_yr{:04d}-{:02d}.png'.format(yr, mo))
   plt.clf()



#plt.show()


