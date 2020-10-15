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
fsz=(14,8)
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-78.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>300.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>310.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #eastern weddell
idx = np.nonzero( (latCell<-65.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<20.0/360.0*2*pii))[0]; size=6;  #entire weddell+
#idx = np.nonzero( (latCell<-60.0/180.0*pii) )[0]; size=1; sz=(14,10); #entire SO
#idx = np.nonzero( (latCell<-50.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<80.0/360.0*2*pii))[0]; size=1.0; fsz=(15,10)  # weddell to Amery 

idxSill=210384-1 # filchner sill

path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'

diffpath = False
diffpath='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist/'

years = np.arange(50,116,1)
months = np.arange(1,13,1)
#months = np.arange(1,2,1)
nt = len(years)*len(months)


fig = plt.figure(1, facecolor='w', figsize=fsz)
nrow=2
ncol=2

diffRange = 0.2

k0=0
k100 = np.argmin(np.absolute(z--100.0))
k50 = np.argmin(np.absolute(z--50.0))
k200 = np.argmin(np.absolute(z--200.0))
k300 = np.argmin(np.absolute(z--300.0))

if diffpath:
   cmap = 'RdBu_r'
else:   
   cmap = 'nipy_spectral'

for yr in years:
 for mo in months:
   print("yr=",yr, "mo=", mo)

   f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   SSS = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,0]
   S50 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k50]
   S100 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k100]
   S200 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k200]
   S300 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k300]

   si=f.variables['timeMonthly_avg_seaIceFreshWaterFlux'][0,idx]
   rain=f.variables['timeMonthly_avg_rainFlux'][0,idx]
   snow=f.variables['timeMonthly_avg_snowFlux'][0,idx]
   evap=f.variables['timeMonthly_avg_evaporationFlux'][0,idx]
   river=f.variables['timeMonthly_avg_riverRunoffFlux'][0,idx]
   iceRunoff=f.variables['timeMonthly_avg_iceRunoffFlux'][0,idx]
   if 'timeMonthly_avg_icebergFreshWaterFlux' in f.variables:
      iceberg=f.variables['timeMonthly_avg_icebergFreshWaterFlux'][0,idx]
   else:
      #print ("Couldn't find timeMonthly_avg_icebergFreshWaterFlux")
      iceberg = si*0.0
   if 'timeMonthly_avg_landIceFreshwaterFlux' in f.variables:
      iceshelf=f.variables['timeMonthly_avg_landIceFreshwaterFlux'][0,idx]
   else:
      #print ("Couldn't find timeMonthly_avg_landIceFreshwaterFlux")
      iceshelf = si*0.0
   FWF=(si+rain+snow+evap+river+iceberg+iceshelf)/1000.0*3.14e7

   u=f.variables['timeMonthly_avg_velocityMeridional'][0,idx,0]
   v=f.variables['timeMonthly_avg_velocityZonal'][0,idx,0]
   #speed=f.variables['timeMonthly_avg_columnIntegratedSpeed'][0,idx]
   speed = (u**2+v**2)*0.5

   f.close()
   
#   f=netCDF4.Dataset('{0}/../../ice/hist/mpascice.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
#   iv=f.variables['timeMonthly_avg_iceVolumeCell'][0,idx]
#   ia=f.variables['timeMonthly_avg_iceAreaCell'][0,idx]
#   f.close()



   if diffpath:
      f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(diffpath, yr, mo), 'r')
      SSSd = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,0]
      Sd50 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k50]
      Sd100 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k100]
      Sd200 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k200]
      Sd300 = f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,k300]
      if 'timeMonthly_avg_landIceFreshwaterFlux' in f.variables:
         iceshelfd=f.variables['timeMonthly_avg_landIceFreshwaterFlux'][0,idx]
      else:
      #print ("Couldn't find timeMonthly_avg_landIceFreshwaterFlux")
         iceshelfd = np.zeros((len(idx),))
      f.close()

   ax0 = fig.add_subplot(nrow, ncol, 1)
   if diffpath:
      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(iceshelf-iceshelfd)/1000.0*3.14e7, vmin=-1.0*1, vmax=1, cmap='RdBu')
   else:
      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(iceshelf), vmin=0, vmax=1, cmap=cmap)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax0)
   ax0.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax0.axis('equal')
#   ax0.set_title("{}:\n year {}, mo {}; ice shelf melt (m/yr)".format(path, yr, mo), fontsize=7)
   ax0.set_title("year {0:04d}, mo {1:02d}; ice shelf melt (m/yr)".format(yr, mo))
   plt.axis('off')


#   ax0 = fig.add_subplot(nrow, ncol, 1)
#   if diffpath:
#      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(iv-ivd), vmin=-1.0*7, vmax=7, cmap='RdBu')
#   else:
#      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(iv), vmin=0, vmax=2, cmap=cmap)
#      sc1.cmap.set_under('k')
#   i2 = np.nonzero(idx==idxSill)[0][0]
#   fig.colorbar(sc1, ax=ax0)
#   ax0.plot(yCell[idxSill], xCell[idxSill], 'k.')
#   ax0.axis('equal')
#   ax0.set_title("{}:\n year {}, mo {}; SI vol (m)".format(path, yr, mo), fontsize=7)
#   plt.axis('off')
#
#
#   ax0 = fig.add_subplot(nrow, ncol, 2)
#   if diffpath:
#      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(FWF-FWFd), vmin=-1.0*7, vmax=7, cmap='RdBu')
#   else:
#      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(FWF), vmin=-10, vmax=10, cmap='RdBu')
#      sc1.cmap.set_under('k')
#   i2 = np.nonzero(idx==idxSill)[0][0]
#   fig.colorbar(sc1, ax=ax0)
#   ax0.plot(yCell[idxSill], xCell[idxSill], 'k.')
#   ax0.axis('equal')
#   #ax0.set_title("{}:\n year {}, mo {}; FWF (m/yr)".format(path, yr, mo), fontsize=7)
#   ax0.set_title("{}:\n year {}, mo {}; FWF (m/yr)".format("", yr, mo), fontsize=7)
#   plt.axis('off')
#
#   ax0 = fig.add_subplot(nrow, ncol, 3)
#   if diffpath:
#      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(speed-speedd), vmin=-1.0*7, vmax=7, cmap=RdBu_r)
#   else:
#      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(speed), vmin=0, vmax=0.01, cmap=cmap)
#      sc1.cmap.set_under('k')
#   i2 = np.nonzero(idx==idxSill)[0][0]
#   fig.colorbar(sc1, ax=ax0)
#   ax0.plot(yCell[idxSill], xCell[idxSill], 'k.')
#   ax0.axis('equal')
#   #ax0.set_title("{}:\n year {}, mo {}; sfc spd (m/s)".format(path, yr, mo), fontsize=7)
#   ax0.set_title("{}:\n year {}, mo {}; sfc spd (m/s)".format("", yr, mo), fontsize=7)
#   plt.axis('off')



   ax0 = fig.add_subplot(nrow, ncol, 2)
   if diffpath:
      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(SSS-SSSd), vmin=-1.0*diffRange, vmax=diffRange, cmap=cmap)
   else:
      sc1 = ax0.scatter(yCell[idx], xCell[idx], s=size, c=(SSS), vmin=33.0, vmax=34.6, cmap=cmap)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax0)
   ax0.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax0.axis('equal')
   #ax0.set_title("{}:\n year {}, mo {}; SSS".format(path, yr, mo), fontsize=7)
   #ax0.set_title("{}:\n year {}, mo {}; SSS".format("", yr, mo), fontsize=7)
   ax0.set_title("year {0:04d}, mo {1:02d}; SSS".format(yr, mo))
   plt.axis('off')



#   ax50 = fig.add_subplot(nrow, ncol, 5)
#   if diffpath:
#      sc1 = ax50.scatter(yCell[idx], xCell[idx], s=size, c=(S50-Sd50), vmin=-1.0*diffRange, vmax=diffRange, cmap=cmap)
#   else:
#      sc1 = ax50.scatter(yCell[idx], xCell[idx], s=size, c=(S50), vmin=33.0, vmax=34.6, cmap=cmap)
#      sc1.cmap.set_under('k')
#   i2 = np.nonzero(idx==idxSill)[0][0]
#   fig.colorbar(sc1, ax=ax50)
#   ax50.plot(yCell[idxSill], xCell[idxSill], 'k.')
#   ax50.axis('equal')
#   #ax50.set_title("{}:\n year {}, mo {}; Sal. 50m".format(path, yr, mo), fontsize=7)
#   ax50.set_title("{}:\n year {}, mo {}; Sal. 50m".format("", yr, mo), fontsize=7)
#   plt.axis('off')

   ax100 = fig.add_subplot(nrow, ncol, 3)
   if diffpath:
      sc1 = ax100.scatter(yCell[idx], xCell[idx], s=size, c=(S100-Sd100), vmin=-1.0*diffRange, vmax=diffRange, cmap=cmap)
   else:
      sc1 = ax100.scatter(yCell[idx], xCell[idx], s=size, c=(S100), vmin=33.0, vmax=34.6, cmap=cmap)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax100)
   ax100.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax100.axis('equal')
   #ax100.set_title("{}:\n year {}, mo {}; Sal. 100m".format(path, yr, mo), fontsize=7)
   ax100.set_title("year {0:04d}, mo {1:02d}; Sal. level {2} (~100m)".format(yr, mo, k100+1))
   plt.axis('off')

   ax200 = fig.add_subplot(nrow, ncol, 4)
   if diffpath:
      sc1 = ax200.scatter(yCell[idx], xCell[idx], s=size, c=(S200-Sd200), vmin=-1.0*diffRange, vmax=diffRange, cmap=cmap)
   else:
      sc1 = ax200.scatter(yCell[idx], xCell[idx], s=size, c=(S200), vmin=33.5, vmax=34.6, cmap=cmap)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=ax200)
   ax200.plot(yCell[idxSill], xCell[idxSill], 'k.')
   ax200.axis('equal')
   #ax200.set_title("{}:\n year {}, mo {}; Sal. 200m".format(path, yr, mo), fontsize=7)
   ax200.set_title("year {0:04d}, mo {1:02d}; Sal. level {2} (~200m)".format(yr, mo, k200))
   plt.axis('off')

#   ax300 = fig.add_subplot(nrow, ncol, 4)
#   if diffpath:
#      sc1 = ax300.scatter(yCell[idx], xCell[idx], s=size, c=(S300-Sd300), vmin=-1.0*diffRange, vmax=diffRange, cmap='RdBu_r')
#   else:
#      sc1 = ax300.scatter(yCell[idx], xCell[idx], s=size, c=(S300), vmin=33.5, vmax=34.6)
#      sc1.cmap.set_under('k')
#   i2 = np.nonzero(idx==idxSill)[0][0]
#   fig.colorbar(sc1, ax=ax300)
#   ax300.plot(yCell[idxSill], xCell[idxSill], 'k.')
#   ax300.axis('equal')
#   ax300.set_title("{}:\n year {}, mo {}; Sal. 300m".format(path, yr, mo), fontsize=7)



   plt.draw()
   plt.savefig('salin_map_yr{:04d}-{:02d}.png'.format(yr, mo))
   plt.clf()



#plt.show()


