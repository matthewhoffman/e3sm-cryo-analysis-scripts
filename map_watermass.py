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
idx = np.nonzero( (latCell<-60.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<80.0/360.0*2*pii))[0]; size=1.4; fsz=(15,6)  # weddell to Amery 

idxSill=210384-1 # filchner sill

path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'

diffpath = False
diffpath='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'

years = np.arange(50,110,1)
months = np.arange(1,13,1)
#months = np.arange(1,2,1)
nt = len(years)*len(months)


fig = plt.figure(1, facecolor='w', figsize=fsz)
nrow=2
ncol=2


def waterMassCalc(f):
   T=f.variables['timeMonthly_avg_activeTracers_temperature']
   S=f.variables['timeMonthly_avg_activeTracers_salinity']
   rho = f.variables['timeMonthly_avg_potentialDensity']
   v  = f.variables['timeMonthly_avg_velocityMeridional']
   
   Tdata=T[0,idx,:]
   Sdata=S[0,idx,:]
   waterMassMask = np.zeros(Tdata.shape)

   waterMassMask[np.where( (Sdata>34.0) * (Tdata>0.0) * (Sdata<40.0) * (Tdata<5.0) )] = 1 # CDW
   waterMassMask[np.where( (Sdata>34.0) * (Tdata<0.0) * (Tdata>-1.5) * (Sdata<40.0 ))] = 2 # mCDW
   waterMassMask[np.where( (Tdata>(-0.0575*Sdata+0.0901)) * (Sdata>34.5) * (Tdata<-1.5) * (Sdata<40.0) )] = 3 # HSSW
   waterMassMask[np.where( (Tdata>(-0.0575*Sdata+0.0901)) * (Sdata<34.5) * (Sdata>34.0) * (Tdata<-1.5) )] = 4 # LSSW
   waterMassMask[np.where( (Tdata>(-0.0575*Sdata+0.0901)) * (Sdata<34.0) * (Sdata>20.0) * (Tdata<5.0)  )] = 5 # AASW
   waterMassMask[np.where( (Tdata<(-0.0575*Sdata+0.0901)) * (Tdata>-3.0) * (Sdata>30.0) * (Sdata<40.0) )] = 6 # ISW

   return waterMassMask


wmLab=('none', 'CDW', 'mCDW', 'HSSW', 'LSSW', 'AASW', 'ISW')

diffRange = 15

for yr in years:
 for mo in months:
   print "yr=",yr, "mo=", mo

   f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   waterMassMask = waterMassCalc(f)
   f.close()

   if diffpath:
      f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(diffpath, yr, mo), 'r')
      waterMassMask2 = waterMassCalc(f)
      f.close()


   axISW = fig.add_subplot(nrow, ncol, 1)
   iMass = 6
   #plt.scatter(yCell[idx], xCell[idx], s=10, c=T[0,idx,:].max(axis=1))
   if diffpath:
      sc1 = axISW.scatter(yCell[idx], xCell[idx], s=size, c=(waterMassMask==iMass).sum(axis=1) - (waterMassMask2==iMass).sum(axis=1), vmin=-1*diffRange, vmax=diffRange, cmap='RdBu')
   else:
      sc1 = axISW.scatter(yCell[idx], xCell[idx], s=size, c=(waterMassMask==iMass).sum(axis=1), vmin=1, vmax=30)
      sc1.cmap.set_under('k')
   i2 = np.nonzero(idx==idxSill)[0][0]
   fig.colorbar(sc1, ax=axISW)
   axISW.plot(yCell[idxSill], xCell[idxSill], 'rx')
   axISW.axis('equal')
   axISW.set_title("{}:\n year {}, mo {}; #layers of {}".format(path, yr, mo, wmLab[iMass]), fontsize=7)

   
   axSW = fig.add_subplot(nrow, ncol, 2)
   if diffpath:
      sc2 = axSW.scatter(yCell[idx], xCell[idx], s=size, c=(np.logical_or(waterMassMask==3, waterMassMask==4)).sum(axis=1) - (np.logical_or(waterMassMask2==3, waterMassMask2==4)).sum(axis=1), vmin=-1*diffRange, vmax=diffRange, cmap='RdBu')
   else:
      sc2 = axSW.scatter(yCell[idx], xCell[idx], s=size, c=(np.logical_or(waterMassMask==3, waterMassMask==4)).sum(axis=1), vmin=1, vmax=60)
      sc2.cmap.set_under([0.5, 0.5, 0.5, 1.0])
   fig.colorbar(sc2, ax=axSW)
   axSW.plot(yCell[idxSill], xCell[idxSill], 'rx')
   axSW.axis('equal')
   axSW.set_title("\nyear {}, mo {}; #layers of SW".format(yr, mo), fontsize=7)

   axAASW = fig.add_subplot(nrow, ncol, 3)
   if diffpath:
      sc3 = axAASW.scatter(yCell[idx], xCell[idx], s=size, c=(waterMassMask==5).sum(axis=1) - (waterMassMask2==5).sum(axis=1), vmin=-1*diffRange, vmax=diffRange, cmap='RdBu')
   else:
      sc3 = axAASW.scatter(yCell[idx], xCell[idx], s=size, c=(waterMassMask==5).sum(axis=1), vmin=1, vmax=30)
      sc3.cmap.set_under([0.5, 0.5, 0.5, 1.0])
   fig.colorbar(sc3, ax=axAASW)
   axAASW.plot(yCell[idxSill], xCell[idxSill], 'rx')
   axAASW.axis('equal')
   axAASW.set_title("\nyear {}, mo {}; #layers of AASW".format(yr, mo), fontsize=7)

   axCDW = fig.add_subplot(nrow, ncol, 4)
   if diffpath:
      sc4 = axCDW.scatter(yCell[idx], xCell[idx], s=size, c=(np.logical_or(waterMassMask==1, waterMassMask==2)).sum(axis=1) - (np.logical_or(waterMassMask2==1, waterMassMask2==2)).sum(axis=1), vmin=-1*diffRange, vmax=diffRange, cmap='RdBu')
   else:
      sc4 = axCDW.scatter(yCell[idx], xCell[idx], s=size, c=(np.logical_or(waterMassMask==1, waterMassMask==2)).sum(axis=1), vmin=1, vmax=60)
      sc4.cmap.set_under([0.5, 0.5, 0.5, 1.0])
   fig.colorbar(sc4, ax=axCDW)
   axCDW.plot(yCell[idxSill], xCell[idxSill], 'rx')
   axCDW.axis('equal')
   axCDW.set_title("\nyear {}, mo {}; #layers of CDW+mCDW".format(yr, mo), fontsize=7)



   plt.draw()
   plt.savefig('water_mass_map_yr{:04d}-{:02d}.png'.format(yr, mo))
   plt.clf()



#plt.show()


