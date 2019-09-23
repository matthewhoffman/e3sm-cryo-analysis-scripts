#!/usr/bin/env python
'''
Script to compare FWF around AIS from different E3SM runs

Prepare climatology data like:
ncra -v timeMonthly_avg_seaIceFreshWaterFlux,timeMonthly_avg_rainFlux,timeMonthly_avg_snowFlux,timeMonthly_avg_evaporationFlux,timeMonthly_avg_riverRunoffFlux,timeMonthly_avg_iceRunoffFlux,timeMonthly_avg_icebergFreshWaterFlux,timeMonthly_avg_landIceFreshwaterFlux,timeMonthly_avg_activeTracers_salinity,timeMonthly_avg_potentialDensity,timeMonthly_avg_salinitySurfaceRestoringTendency mpaso.hist.am.timeSeriesStatsMonthly.002* /global/cscratch1/sd/hoffman2/acme_scratch/analysis/mpaso.hist.am.timeSeriesStatsMonthly.mean.0020-0029.G-IAF-ISMF.nc
'''

from __future__ import absolute_import, division, print_function, unicode_literals

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
areaCell = fmesh.variables['areaCell'][:]
z = np.zeros(depths.shape)
z[0] = -0.5 * depths[0]
z[1:] = -0.5 * (depths[0:-1] + depths[1:])
pii=3.14159

# ---- Choose spatial extent -----
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-77.2/180.0*pii, latCell>-78.5/180.0*pii), np.logical_and(lonCell>317.0/360.0*2.0*pii, lonCell<326.0/360.0*2*pii)))[0]  #Darelius 2016 fig 3
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-78.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>300.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>270.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell wider
#idx = np.nonzero( (latCell<-60.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<80.0/360.0*2*pii))[0]; size=1.4; fsz=(15,9)  # weddell to Amery 
idx = np.nonzero( (latCell<-50.0/180.0*pii) )[0] # SO
size = 0.5
# -------------------------



dir='/global/cscratch1/sd/hoffman2/acme_scratch/analysis/'
flist=(
        'mpaso.hist.am.timeSeriesStatsMonthly.mean.0080-0089.G-IAF.nc',
        'mpaso.hist.am.timeSeriesStatsMonthly.mean.0020-0029.G-IAF-DIB.nc',
        'mpaso.hist.am.timeSeriesStatsMonthly.mean.0081-0089.G-IAF-ISMF.nc',
        'mpaso.hist.am.timeSeriesStatsMonthly.mean.0020-0029.G-IAF-DIB-ISMF.nc'
        , 'mpaso.hist.am.timeSeriesStatsMonthly.mean.0020-0029.B-DIB.nc'
        , 'mpaso.hist.am.timeSeriesStatsMonthly.mean.0020-0029.B-DIB-ISMF.nc'
        )


fig1 = plt.figure(1, facecolor='w') # TS
nrow=len(flist)
ncol=2

fn=0
for fname in flist:
   print("\n----------")
   print(fname)
   f=netCDF4.Dataset(dir+fname, 'r')
   si=f.variables['timeMonthly_avg_seaIceFreshWaterFlux'][0,idx]
   rain=f.variables['timeMonthly_avg_rainFlux'][0,idx]
   snow=f.variables['timeMonthly_avg_snowFlux'][0,idx]
   evap=f.variables['timeMonthly_avg_evaporationFlux'][0,idx]
   river=f.variables['timeMonthly_avg_riverRunoffFlux'][0,idx]
   iceRunoff=f.variables['timeMonthly_avg_iceRunoffFlux'][0,idx]
   if 'timeMonthly_avg_icebergFreshWaterFlux' in f.variables:
      iceberg=f.variables['timeMonthly_avg_icebergFreshWaterFlux'][0,idx]
   else:
      print ("Couldn't find timeMonthly_avg_icebergFreshWaterFlux")
      iceberg = si*0.0
   if 'timeMonthly_avg_landIceFreshwaterFlux' in f.variables:
      iceshelf=f.variables['timeMonthly_avg_landIceFreshwaterFlux'][0,idx]
   else:
      print ("Couldn't find timeMonthly_avg_landIceFreshwaterFlux")
      iceshelf = si*0.0
#   sal=f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,0]
#   dens=f.variables['timeMonthly_avg_potentialDensity'][0,idx,0]
#   restTend=f.variables['timeMonthly_avg_salinitySurfaceRestoringTendency'][0,idx]


   AISFWF = river + iceRunoff + iceberg + iceshelf
   if fn==0:
     refAISFWF=AISFWF

   if fn==0:
       ax1 = fig1.add_subplot(nrow, ncol, fn*2+1); ax1.set_aspect('equal')
   else:
       ax = fig1.add_subplot(nrow, ncol, fn*2+1, sharex=ax1, sharey=ax1); ax.set_aspect('equal')
   plt.title(fname); plt.axis('off')
   rng = max(abs(AISFWF.min()), abs(AISFWF.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=AISFWF, vmin=0.0, vmax=0.00002)
   plt.colorbar()
   
   ax = fig1.add_subplot(nrow, ncol, fn*2+2, sharex=ax1, sharey=ax1); ax.set_aspect('equal')
   plt.title('diff'); plt.axis('off')
#   rng = max(abs(AISFWF.min()), abs(AISFWF.max()))
#   vmin = -1.0*rng; vmax=rng
   rng = 0.000004
   plt.scatter(yCell[idx], xCell[idx], s=size, c=AISFWF-refAISFWF, vmin=-1*rng, vmax=rng, cmap='RdBu')
   plt.colorbar()

   plt.draw()

   # global stats
   areaSum=areaCell[idx].sum()
#   factor = 1.0/1000.0*3.14e7; print('regional budget terms (m/yr)')
#   factor = areaSum/1.0e6/1000.0; print('regional budget terms (fwSv)')
   factor = areaSum/1.0e12*3.14e7; print('* regional budget terms (Gt/yr):')
   print("iceberg={}".format((iceberg*areaCell[idx]).sum()/areaSum*factor))
   print("iceshelf={}".format((iceshelf*areaCell[idx]).sum()/areaSum*factor))
   print("river={}".format((river*areaCell[idx]).sum()/areaSum*factor))
   print("iceRunoff={}".format((iceRunoff*areaCell[idx]).sum()/areaSum*factor))
   print("total={}".format((AISFWF*areaCell[idx]).sum()/areaSum*factor))


   f.close()

   fn+=1

plt.show()
