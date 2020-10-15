#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
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

fmesh=netCDF4.Dataset('/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
areaCell = fmesh.variables['areaCell'][:]


pii=3.14159


# ---- Choose time(s) -----
yrs=(95,)
#yrs = np.arange(95,102,1)
mos=(1,)
#mos=np.arange(1,13,1)
# -------------------------


# ---- Choose spatial extent -----
idx = np.nonzero(np.logical_and(np.logical_and(latCell<-77.4/180.0*pii, latCell>-78.2/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<327.0/360.0*2*pii)))[0]  #Darelius 2016 fig 3
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-78.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>300.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>270.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell wider
#idx = np.nonzero( (latCell<-60.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<80.0/360.0*2*pii))[0]; size=1.4; fsz=(15,9)  # weddell to Amery 
#idx = np.nonzero( (latCell<-50.0/180.0*pii) )[0] # SO
# -------------------------


# ---- Choose run directory ----
path='/project/projectdirs/m3412/simulations/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'
#path='/project/projectdirs/m3412/simulations/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist/'
# -------------------------


print("Found {} cells".format(len(idx)))

figMap = plt.figure(100, facecolor='w')
idx2 = np.nonzero(latCell<-60.0/180.0*pii)[0]
plt.plot(yCell[idx2], xCell[idx2], 'k.')
plt.plot(yCell[idx], xCell[idx], 'r.')


fig1 = plt.figure(1, facecolor='w') # TS
nrow=1
ncol=1
doTS=True

fig2 = plt.figure(2, facecolor='w') # velo
doVelo=False

figFW = plt.figure(3, facecolor='w') # FW sfc budget
#figFW2 = plt.figure(4, facecolor='w') # FW sfc budget
doFW=True
doFW=False

nVertLev = len(fmesh.dimensions['nVertLevels'])
zs = np.zeros((len(idx), nVertLev))
cnt=0
for i in idx:
   maxLevelCell=fmesh.variables['maxLevelCell'][i]
   bottomDepth = fmesh.variables['bottomDepth'][i]
   layerThickness = fmesh.variables['layerThickness'][0, i, :maxLevelCell]
   thicknessSum = layerThickness.sum()
   thicknessCumSum = layerThickness.cumsum()
   zSurf = bottomDepth - thicknessSum
   zLayerBot = zSurf - thicknessCumSum
   z = zLayerBot + 0.5 * layerThickness
   zs[cnt, :maxLevelCell] = z
   cnt += 1



for yr in yrs:
 for mo in mos:
  print("yr=",yr, "mo=", mo)
  f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
#  f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.mean.0050-0059.nc'.format(path), 'r')
  #f=netCDF4.Dataset('/global/cscratch1/sd/hoffman2/acme_scratch/analysis/mpaso.hist.am.timeSeriesStatsMonthly.mean.0080-0089.G-IAF.nc', 'r')

  if doTS:
   fig1.clf()
   axTS = fig1.add_subplot(nrow, ncol, 1)
   plt.sca(axTS)
   plt.ylabel('temperature (deg. C)')
   plt.xlabel('salinity (psu)')
   plt.grid()


   T=f.variables['timeMonthly_avg_activeTracers_temperature']
   S=f.variables['timeMonthly_avg_activeTracers_salinity']
   
   # get data
   Ts = T[0,idx,:]
   Ss = S[0,idx,:]
#   Ts[zs>-200.0]=0.0; Ss[zs>-200.0]=0.0 # remove shallow values
   cmap = plt.cm.get_cmap('viridis', 64)
   cmap.set_over('0.8')
   sc=axTS.scatter(Ss[:], Ts[:], s=1, c=zs, vmin=-700, vmax=-200.0, cmap=cmap)
   
   #axTS.set_ylim([-2.6,2.0])
   #axTS.set_xlim([33.0,34.8])
   axTS.set_ylim([-2.35,-1.4]); axTS.set_xlim([34.1,34.75])  # Darelius 2016 range
   plt.colorbar(sc)
   
   plt.plot([34.0, 34.0], [-1.85, 0.2], 'k:')
   plt.plot([34.5, 34.5], [-1.86, -1.5], 'k:')
   plt.plot([34.0, 35.0], [-1.5, -1.5], 'k:')
   plt.plot([34.0, 35.0], [0.0, 0.0], 'k:')
   plt.plot([33.5, 35.0], [-0.0575*33.5+0.0901, -0.0575*35.0+0.0901], 'k:')
   plt.title("year {}, mo {}".format(yr, mo))
   plt.draw()
   plt.savefig('TS_yr{:04d}-{:02d}.png'.format(yr, mo))
   
   
   
   # velo plots ------------------
  if doVelo:
   v=f.variables['timeMonthly_avg_velocityMeridional'][0,idx,:]
   u=f.variables['timeMonthly_avg_velocityZonal'][0,idx,:]
   swsv=f.variables['timeMonthly_avg_windStressMeridional'][0,idx]
   swsu=f.variables['timeMonthly_avg_windStressZonal'][0,idx]

   print(" calc mean velo")
   maxLevelCell=fmesh.variables['maxLevelCell'][idx]
   layerThickness = fmesh.variables['layerThickness'][0, idx, :]
   umean=np.zeros((len(idx),))
   vmean=np.zeros((len(idx),))
   for i in range(len(idx)):
       layersum=0.0
       maxLevelHere = maxLevelCell[i]
       umean[i]=(u[i,:maxLevelHere]*layerThickness[i,:maxLevelHere]).sum()/layerThickness[i,:maxLevelHere].sum()
       vmean[i]=(v[i,:maxLevelHere]*layerThickness[i,:maxLevelHere]).sum()/layerThickness[i,:maxLevelHere].sum()
   print("done")
   umag = (umean**2+vmean**2)**0.5
 
   fig2.clf()
   
   axUsrf = fig2.add_subplot(2, 2, 1)
   axUsrf.axis('equal')
   plt.quiver(yCell[idx], xCell[idx], u[:,0], v[:,0])
   print("d={}, z={}".format(0, z[0]))
   plt.title("year={}, depth={}m".format(yr, z[0]))
   plt.title
   
   
   d=np.argmin(np.absolute(z--250.0))
   print("d={}, z={}".format(d, z[d]))
   
   axU2 = fig2.add_subplot(2, 2, 2)
   axU2.axis('equal')
   plt.quiver(yCell[idx], xCell[idx], u[:,d], v[:,d])
   plt.title("year={}, depth={}m".format(yr, z[d]))

   axUmn = fig2.add_subplot(2,2, 3)
   axUmn.axis('equal')
   #plt.tricontourf(yCell[idx], xCell[idx], umag); plt.colorbar()
   plt.quiver(yCell[idx], xCell[idx], umean[:]/umag, vmean[:]/umag,  np.minimum(umag, 0.03),  units='width')
   #plt.streamplot(Yi, Xi, ureg, vreg)#, density=[0.5, 1])
   plt.title("year={}, mean velo".format(yr))

 
   axWS = fig2.add_subplot(2, 2, 4)
   axWS.axis('equal')
   plt.quiver(yCell[idx], xCell[idx], swsu[:], swsv[:])
   plt.title("year={}, surface stress".format(yr))
   plt.draw()
   plt.savefig('velo_yr{:04d}.png'.format(yr))


   # ------- FW surface budget ----
  if doFW:
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
   sal=f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,0]
   dens=f.variables['timeMonthly_avg_potentialDensity'][0,idx,0]
   restTend=f.variables['timeMonthly_avg_salinitySurfaceRestoringTendency'][0,idx]

   PE = rain + snow + evap
   restoring = -1.0*restTend/sal*dens
   total = PE + si + iceRunoff + iceshelf + restoring + iceberg + river

   
   figFW.clf()
   nrow = 4
   ncol = 2
   size = 0.5
   vmin = min(PE.min(), si.min(), iceRunoff.min(), iceshelf.min(), restoring.min(), river.min())
   vmax = max(PE.max(), si.max(), iceRunoff.max(), iceshelf.max(), restoring.max(), river.max())

   axPE = figFW.add_subplot(nrow, ncol, 1); axPE.axis('equal')
   plt.title('P+S-E'); plt.axis('off')
   rng = max(abs(PE.min()), abs(PE.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=PE, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()

   axSI = figFW.add_subplot(nrow, ncol, 2); axSI.axis('equal')
   plt.title('sea ice'); plt.axis('off')
   rng = max(abs(si.min()), abs(si.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=si, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()

   axRun = figFW.add_subplot(nrow, ncol, 3); axRun.axis('equal')
   plt.title('iceRunoff'); plt.axis('off')
   rng = max(abs(iceRunoff.min()), abs(iceRunoff.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=iceRunoff, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()

   axShelf = figFW.add_subplot(nrow, ncol, 4); axShelf.axis('equal')
   plt.title('iceshelf basal melt'); plt.axis('off')
   rng = max(abs(iceshelf.min()), abs(iceshelf.max()))
   vmin = -1.0*rng; vmax=rng
   vmin = -0.00004; vmax=0.00004
   plt.scatter(yCell[idx], xCell[idx], s=size, c=iceshelf, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()

   axRestor = figFW.add_subplot(nrow, ncol, 5); axRestor.axis('equal')
   plt.title('restoring'); plt.axis('off')
   rng = max(abs(restoring.min()), abs(restoring.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=restoring, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()

   axIB = figFW.add_subplot(nrow, ncol, 6); axIB.axis('equal')
   plt.title('iceberg'); plt.axis('off')
   rng = max(abs(iceberg.min()), abs(iceberg.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=iceberg, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()
 
   axRiver = figFW.add_subplot(nrow, ncol, 7); axRiver.axis('equal')
   plt.title('river'); plt.axis('off')
   rng = max(abs(river.min()), abs(river.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=river, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()

   axFW = figFW.add_subplot(nrow,ncol,8); axFW.axis('equal')
   plt.title('total FW'); plt.axis('off')
   rng = max(abs(total.min()), abs(total.max()))
   vmin = -1.0*rng; vmax=rng
   plt.scatter(yCell[idx], xCell[idx], s=size, c=total, vmin=vmin, vmax=vmax, cmap='RdBu')
   plt.colorbar()

   plt.draw()
 
   # global stats
   areaSum=areaCell[idx].sum()
   factor = 1.0/1000.0*3.14e7; print('regional budget terms (m/yr)')
   factor = areaSum/1.0e6/1000.0; print('regional budget terms (fwSv)')
   print("P+S={}".format(((rain+snow)*areaCell[idx]).sum()/areaSum*factor))
   print("E={}".format((evap*areaCell[idx]).sum()/areaSum*factor))
   print("PE={}".format((PE*areaCell[idx]).sum()/areaSum*factor))
   print("sea ice={}".format((si*areaCell[idx]).sum()/areaSum*factor))
   print("restoring={}".format((restoring*areaCell[idx]).sum()/areaSum*factor))
   print("iceberg={}".format((iceberg*areaCell[idx]).sum()/areaSum*factor))
   print("iceshelf={}".format((iceshelf*areaCell[idx]).sum()/areaSum*factor))
   print("river={}".format((river*areaCell[idx]).sum()/areaSum*factor))
   print("total={}".format((total*areaCell[idx]).sum()/areaSum*factor))


if len(yrs)==1 and len(mos)==1:
   #"interative" mode
   plt.show()
