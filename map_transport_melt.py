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
import matplotlib

fmesh=netCDF4.Dataset('/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
areaCell = fmesh.variables['areaCell'][:]


pii=3.14159


# ---- Choose time(s) -----
yrs=(50,)
#yrs = np.arange(95,102,1)
mos=(7,)
#mos=np.arange(1,13,1)
# -------------------------


# ---- Choose spatial extent -----
idx = np.nonzero(np.logical_and(np.logical_and(latCell<-77.4/180.0*pii, latCell>-78.2/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<327.0/360.0*2*pii)))[0]  #Darelius 2016 fig 3
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-75.0/180.0*pii, latCell>-78.0/180.0*pii), np.logical_and(lonCell>320.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner trough in front of ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-78.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>315.0/360.0*2.0*pii, lonCell<330.0/360.0*2*pii)))[0]  #filchner ice shelf
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>300.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell
idx = np.nonzero(np.logical_and(np.logical_and(latCell<-70.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>270.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  #entire weddell wider
#idx = np.nonzero( (latCell<-60.0/180.0*pii) * (latCell>-85.0/180.0*pii) * np.logical_or(lonCell>280.0/360.0*2.0*pii, lonCell<80.0/360.0*2*pii))[0]; size=1.4; fsz=(15,9)  # weddell to Amery 
#idx = np.nonzero(np.logical_and(np.logical_and(latCell<-72.0/180.0*pii, latCell>-85.0/180.0*pii), np.logical_and(lonCell>291.0/360.0*2.0*pii, lonCell<350.0/360.0*2*pii)))[0]  # Daae 2020 region
#idx = np.nonzero( (latCell<-50.0/180.0*pii) )[0] # SO
# -------------------------


# ---- Choose run directory ----
path='/project/projectdirs/m3412/simulations/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'
#path='/project/projectdirs/m3412/simulations/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist/'
# -------------------------


print("Found {} cells".format(len(idx)))

#figMap = plt.figure(100, facecolor='w')
#idx2 = np.nonzero(latCell<-60.0/180.0*pii)[0]
#plt.plot(yCell[idx2], xCell[idx2], 'k.')
#plt.plot(yCell[idx], xCell[idx], 'r.')


figTransportMelt = plt.figure(4, facecolor='w') # velo

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

  v=f.variables['timeMonthly_avg_velocityMeridional'][0,idx,:]
  u=f.variables['timeMonthly_avg_velocityZonal'][0,idx,:]
  swsv=f.variables['timeMonthly_avg_windStressMeridional'][0,idx]
  swsu=f.variables['timeMonthly_avg_windStressZonal'][0,idx]

  melt = f.variables['timeMonthly_avg_landIceFreshwaterFlux'][0,idx] / 910.0 * 3600.0 * 24.0 * 365.0

  print(" calc mean velo")
  maxLevelCell=fmesh.variables['maxLevelCell'][idx]
  layerThickness = fmesh.variables['layerThickness'][0, idx, :]
  umean=np.zeros((len(idx),))
  vmean=np.zeros((len(idx),))
  uBtm=np.zeros((len(idx),))
  vBtm=np.zeros((len(idx),))
  for i in range(len(idx)):
      layersum=0.0
      maxLevelHere = maxLevelCell[i]
      umean[i]=(u[i,:maxLevelHere]*layerThickness[i,:maxLevelHere]).sum()/layerThickness[i,:maxLevelHere].sum()
      vmean[i]=(v[i,:maxLevelHere]*layerThickness[i,:maxLevelHere]).sum()/layerThickness[i,:maxLevelHere].sum()
      uBtm[i]=u[i,maxLevelHere-1]
      vBtm[i]=v[i,maxLevelHere-1]
  print("done")
  umag = (umean**2+vmean**2)**0.5
  umagBtm = (uBtm**2+vBtm**2)**0.5
 

  # Transport and melt plot
  axUmn = figTransportMelt.add_subplot(1,1,1)
  axUmn.axis('equal')
  #plt.tricontourf(yCell[idx], xCell[idx], umag); plt.colorbar()
  umean[umag<0.005]=0.0; vmean[umag<0.005]=0.0
  #plt.tricontourf(yCell[idx], xCell[idx], melt,  levels=np.linspace(-1.5, 1.5, 10), cmap='bwr'); plt.colorbar()
  plt.scatter(yCell[idx], xCell[idx], s=50, c=melt,  vmin=-2, vmax=2, cmap='bwr'); plt.colorbar() # bwr
  plt.quiver(yCell[idx], xCell[idx], umean[:]/umag, vmean[:]/umag, np.minimum(1.0e-1,np.maximum(0.5e-2,umag)), norm=matplotlib.colors.LogNorm(), cmap='Greys'); plt.colorbar()  #Greys
#  plt.quiver(yCell[idx], xCell[idx], umean[:], vmean[:],  np.maximum(-2, np.log10(umag)), cmap='Greys'); plt.colorbar()
  #plt.streamplot(Yi, Xi, ureg, vreg)#, density=[0.5, 1])
  plt.title("year={}, mean velo".format(yr))

 
  plt.draw()

if len(yrs)==1 and len(mos)==1:
   #"interative" mode
   plt.show()
