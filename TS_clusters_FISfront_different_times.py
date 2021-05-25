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
import gsw
from gsw.density import sigma0

fmesh=netCDF4.Dataset('/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
areaCell = fmesh.variables['areaCell'][:]
bottomDepth = fmesh.variables['bottomDepth'][:]


pii=3.14159


# ---- Choose time(s) -----
#yrs=(10, 30, 50, 60, 70,80, 90, 100, 150)
#yrs=(10, 30, 50, 60, 70,80, 90, 100, 100)
#yrs=np.arange(1,150,1)
yrs=np.arange(1,120,2)

#yrs=(50, 50, 50, 60, 70,80, 90, 100, 150)
#yrs = np.arange(95,102,1)
mos=(1,)

mo=1
#mos=np.arange(1,13,1)
# -------------------------

# ---- Choose run directory ----
path='/project/projectdirs/m3412/simulations/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'
path='/project/projectdirs/m3412/simulations/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist/'
path='/project/projectdirs/m3412/simulations/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/archive/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/e3sm_scratch/cori-knl/20191003.GMPAS-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist'
# -------------------------



# ---- Choose spatial extent -----
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-77.4/180.0*pii, latCell>-78.2/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<327.0/360.0*2*pii)))[0]  #Darelius 2016 fig 3
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-74.0/180.0*pii, latCell>-76.0/180.0*pii), np.logical_and(lonCell>(360.0-65.0)/360.0*2.0*pii, lonCell<(360.0-55.0)/360.0*2*pii)))[0]  # Ronne inflow
idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-78.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
idx = np.nonzero( (latCell<-77.5/180.0*pii) * (latCell>-78.0/180.0*pii) * (lonCell>320.0/360.0*2.0*pii) * (lonCell<330.0/360.0*2*pii) * (bottomDepth>500.0))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-76.0/180.0*pii), np.logical_and(lonCell>325.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>300.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>270.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell wider
#idx = np.nonzero( (latCell<-60.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<80.0/360.0*2*pii))[0]; size=1.4; fsz=(15,9)  # weddell to Amery 
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-72.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>291.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  # Daae 2020 region
#idx = np.nonzero( (latCell<-50.0/180.0*pii) )[0] # SO

# -------------------------

idx2 = np.nonzero(np.logical_and(np.logical_and(latCell<-73.0/180.0*pii, latCell>-74.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<335.0/360.0*2*pii)))[0]  # offshore of Filcher Trough
idx2 = np.nonzero(np.logical_and(np.logical_and(latCell<-74.0/180.0*pii, latCell>-74.5/180.0*pii), np.logical_and(lonCell>325.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  # offshore of Filcher Trough - tighter
idx2 = np.nonzero( (latCell<-72/180.0*pii) * (latCell>-75.1/180.0*pii) * (lonCell>325.0/360.0*2.0*pii) * (lonCell<335.0/360.0*2*pii) * (bottomDepth>700.0) * (bottomDepth<2500))[0]  #filchner trough in front of ice shelf
# -------------------------

print("Found {} cells in idx".format(len(idx)))
print("Found {} cells in idx2".format(len(idx2)))

figMap = plt.figure(100, facecolor='w')
idxSO = np.nonzero(latCell<-60.0/180.0*pii)[0]
plt.plot(yCell[idxSO], xCell[idxSO], 'k.')
plt.plot(yCell[idx], xCell[idx], 'r.')
plt.plot(yCell[idx2], xCell[idx2], 'b.')


fig1 = plt.figure(1, facecolor='w') # TS
nrow=1
ncol=1
axTS = fig1.add_subplot(nrow, ncol, 1)
plt.sca(axTS)
plt.ylabel('temperature (deg. C)')
plt.xlabel('salinity (psu)')
#plt.grid()

c='0.85'; lw=2
plt.plot([34.0, 34.0], [-1.85, 0.2], '-', color=c, linewidth=lw, zorder=1)
plt.plot([34.5, 34.5], [-1.86, -1.5], '-', color=c, linewidth=lw, zorder=1)
plt.plot([34.0, 35.0], [-1.5, -1.5], '-', color=c, linewidth=lw, zorder=1)
plt.plot([34.0, 35.0], [0.0, 0.0], '-', color=c, linewidth=lw, zorder=1)
plt.plot([33.5, 35.0], [-0.0575*33.5+0.0901, -0.0575*35.0+0.0901], '-', color=c, linewidth=lw, zorder=1)

nVertLev = len(fmesh.dimensions['nVertLevels'])
zs = np.zeros((len(idx), nVertLev))
cnt=0
for i in idx:
   maxLevelCell=fmesh.variables['maxLevelCell'][i]
   bottomDepthHere = bottomDepth[i]
   layerThickness = fmesh.variables['layerThickness'][0, i, :maxLevelCell]
   thicknessSum = layerThickness.sum()
   thicknessCumSum = layerThickness.cumsum()
   zSurf = bottomDepthHere - thicknessSum
   zLayerBot = zSurf - thicknessCumSum
   z = zLayerBot + 0.5 * layerThickness
   zs[cnt, :maxLevelCell] = z
   cnt += 1

zs2 = np.zeros((len(idx2), nVertLev))
cnt=0
for i in idx2:
   maxLevelCell=fmesh.variables['maxLevelCell'][i]
   bottomDepthHere = bottomDepth[i]
   layerThickness = fmesh.variables['layerThickness'][0, i, :maxLevelCell]
   thicknessSum = layerThickness.sum()
   thicknessCumSum = layerThickness.cumsum()
   zSurf = bottomDepthHere - thicknessSum
   zLayerBot = zSurf - thicknessCumSum
   z = zLayerBot + 0.5 * layerThickness
   zs2[cnt, :maxLevelCell] = z
   cnt += 1


# Build arrays for each cluster
cluster1 = np.zeros((2, len(yrs)))
cluster2 = np.zeros((2, len(yrs)))
i=0
for yr in yrs:
# for mo in mos:
  print("yr=",yr)#, "mo=", mo)
  try:
    f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
  except:
    continue

  # cluster1
  Ts=f.variables['timeMonthly_avg_activeTracers_temperature'][0,idx,:]
  Ss=f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,:]
  badrng=(zs>-200.0)
  Ts[badrng]=0.0; Ss[badrng]=0.0 # remove shallow values
  cluster1[:,i] = [Ss[Ss>30].flatten().mean(), Ts[Ts!=0].flatten().mean()]
  # cluster2
  Ts=f.variables['timeMonthly_avg_activeTracers_temperature'][0,idx2,:]
  Ss=f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx2,:]
  badrng=np.logical_or(zs2>-200.0, zs2<-800.0)
  Ts[badrng]=0.0; Ss[badrng]=0.0 # remove shallow & deep values
  cluster2[:,i] = [Ss[Ss>30].flatten().mean(), Ts[Ts!=0].flatten().mean()]
  i += 1

#Plot
#density contours
Tbins = np.arange(-5,5,0.1)
Sbins = np.arange(33.8, 35, 0.02)
PT, SP = np.meshgrid(Tbins, Sbins)
SA = gsw.SA_from_SP(SP, p=0., lon=0., lat=-75.)
CT = gsw.CT_from_t(SA, PT, p=0.)
neutralDensity = sigma0(SA, CT)
rhoInterval = 0.05
contours = np.arange(24., 29.+rhoInterval, rhoInterval)
CS = plt.contour(SP, PT, neutralDensity, contours, linewidths=1.,
                 colors='0.8', zorder=2)
#plt.clabel(CS, fontsize=12, inline=1, fmt='%4.2f')

colors=[ cm.jet(x) for x in np.linspace(0.0, 1.0, len(yrs)) ] 
cmap = plt.cm.get_cmap('jet')
#colors=[ cm.jet(x) for x in np.linspace(0.0, 1.0, 152) ] 
#  plt.plot(Ss[:].flatten(), Ts[:].flatten(), '.', markersize=2, color=colors[iCol], label=yr)
plt.plot(cluster1[0,:], cluster1[1,:], '--', color=[0.5, 0.5, 0.5], zorder=5)
sc=plt.scatter(cluster1[0,:], cluster1[1,:], s=20, c=yrs, zorder=10, cmap=cmap)
plt.colorbar(sc)
#plt.scatter(cluster1[0,:], cluster1[1,:], c=colors,  markersize=4)
if yr in np.arange(63,68) or yr in np.arange(63+62, 63+62+5):
   plt.plot(cluster1[0,yr], cluster2[1,yr], 'x', markersize=4, color='k', label=yr, zorder=8)
  
plt.scatter(cluster2[0,:], cluster2[1,:], s=10, c='w', edgecolors=colors, marker='^', zorder=10, cmap=cmap)
#plt.plot(cluster2[0,:], cluster2[1,:], 'xk', markersize=2)
'''
path='/project/projectdirs/m3412/simulations/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/archive/ocn/hist'
colors=[ cm.jet(x) for x in np.linspace(0.0, 1.0, len(yrs)) ] 
iCol=0
for yr in yrs:
# for mo in mos:
  print("yr=",yr)#, "mo=", mo)
  f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')

  Ts=f.variables['timeMonthly_avg_activeTracers_temperature'][0,idx,:]
  Ss=f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,:]
  Ts[zs>-200.0]=0.0; Ss[zs>-200.0]=0.0 # remove shallow values
#  plt.plot(Ss[:].flatten(), Ts[:].flatten(), '.', markersize=2, color=colors[iCol], label=yr)
  plt.plot(Ss[Ss>30].flatten().mean(), Ts[Ts!=0].flatten().mean(), '^', markersize=10, color=colors[iCol])
  iCol+=1
  
path='/project/projectdirs/m3412/simulations/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist/'
yrs=(50, 50, 50, 60, 70,80, 90, 100, 150)
colors=[ cm.jet(x) for x in np.linspace(0.0, 1.0, len(yrs)) ] 
iCol=0
for yr in yrs:
# for mo in mos:
  print("yr=",yr)#, "mo=", mo)
  f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')

  Ts=f.variables['timeMonthly_avg_activeTracers_temperature'][0,idx,:]
  Ss=f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,:]
  Ts[zs>-200.0]=0.0; Ss[zs>-200.0]=0.0 # remove shallow values
#  plt.plot(Ss[:].flatten(), Ts[:].flatten(), '.', markersize=2, color=colors[iCol], label=yr)
  plt.plot(Ss[Ss>30].flatten().mean(), Ts[Ts!=0].flatten().mean(), 'o', markersize=10, color=colors[iCol])
  
  iCol+=1
''' 

axTS.set_ylim([-2.2,1.7]); axTS.set_xlim([34.1,34.8])
#axTS.set_ylim([-2.35,-1.4]); axTS.set_xlim([34.1,34.75])  # Darelius 2016 range
#axTS.set_ylim([-3,1.5]); axTS.set_xlim([34.0,34.8])  # Daae 2020 S6 range

#plt.colorbar(sc)

#axTS.legend()
#plt.colorbar()
plt.draw()
plt.savefig('TS_cluster.png')
   
   
plt.show()
