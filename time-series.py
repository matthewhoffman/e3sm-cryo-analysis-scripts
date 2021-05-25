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

fmesh=netCDF4.Dataset('/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
areaCell= fmesh.variables['areaCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
z = np.zeros(depths.shape)
z[0] = -0.5 * depths[0]
z[1:] = -0.5 * (depths[0:-1] + depths[1:])
pii=3.14159



idxBIS = np.nonzero( (latCell<-71.0/180.0*pii) * (latCell>-76.0/180.0*pii) * (lonCell>330.0/360.0*2.0*pii) * (lonCell<350.0/360.0*2*pii) )[0]
idxFS =  np.nonzero( (latCell<-74.0/180.0*pii) * (latCell>-76.0/180.0*pii) * (lonCell>320.0/360.0*2.0*pii) * (lonCell<330.0/360.0*2*pii) )[0]
idxFT =  np.nonzero( (latCell<-76.0/180.0*pii) * (latCell>-78.0/180.0*pii) * (lonCell>320.0/360.0*2.0*pii) * (lonCell<330.0/360.0*2*pii) )[0]


idxEAIS = np.nonzero(  (lonCell>330.0/360.0*2.0*pii) + (lonCell<60.0/360.0*2*pii) )[0]
#idxSill=210384-1 # filchner sill

path='/project/projectdirs/m3412/simulations/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'

#diffpath = False
#diffpath='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'

#years = np.arange(50,110,1)
years = np.arange(50,55,1)
months = np.arange(1,13,1)
months = np.arange(1,2,1)
nt = len(years)*len(months)




k0 = np.nonzero( (z<0.0) * (z>-100.0))[0]
k100 = np.nonzero( (z<-100.0) * (z>-200.0))[0]
k200 = np.nonzero( (z<-200.0) * (z>-300.0))[0]
k300 = np.nonzero( (z<-300.0) * (z>-400.0))[0]
k400 = np.nonzero( (z<-400.0) * (z>-500.0))[0]
k500 = np.nonzero( (z<-500.0) )[0]

Zbins = (k0, k100, k200, k300, k400)
nZbins = len(Zbins)
ZbinName=('0-100', '100-200', '200-300', '300-400', '400-500', '500+')

Boxes = (idxBIS, idxFS, idxFT)
nBoxes = len(Boxes)
BoxName = ('Brunt IS', 'Filchner Sill', 'Filchner Trough')

Sall = np.zeros((nt, nBoxes, nZbins))
Tall = np.zeros((nt, nBoxes, nZbins))

Mall = np.zeros((nt, nBoxes))
MEAIS = np.zeros((nt,))

times = np.zeros((nt,))


t=0
for yr in years:
 for mo in months:
   print("yr=",yr, "mo=", mo)

   times[t] = yr+(mo-1.0)/12.0

   f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   Mdata = f.variables['timeMonthly_avg_landIceFreshwaterFlux'][0, idxEAIS]
   MEAIS[t] = (Mdata*areaCell[idxEAIS]).sum()
   for i in range(nBoxes):
      Mdata = f.variables['timeMonthly_avg_landIceFreshwaterFlux'][0, Boxes[i]]
      Mall[t,i] = (Mdata*areaCell[Boxes[i]]).sum()
      for k in range(nZbins):
         Sdata = f.variables['timeMonthly_avg_activeTracers_salinity'][0, Boxes[i], Zbins[k]]
         Sall[t, i, k] = np.ma.array(Sdata, mask=( (Sdata<30.0) + (Sdata>40.0) ) ).mean()
   f.close()

   t += 1


print("beginning plot")
fsz=(12,9)
fig = plt.figure(1, facecolor='w', figsize=fsz)
nrow=nBoxes+1
ncol=1

axM = fig.add_subplot(nrow, ncol, 1)
plt.plot(times, Mall[:,0]*-5, 'k-', label="Brunt x5")
plt.plot(times, MEAIS[:]*-1, 'b-', label="EAIS")
plt.ylabel('relative melt (flipped)')
plt.title('melt')
plt.legend()

colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, nZbins) ]
styles = ('solid','dashed', 'dotted')

for i in range(nBoxes):
   ax = fig.add_subplot(nrow, ncol, i+2, sharex=axM)
   for k in range(nZbins):
      plt.plot(times, Sall[:,i,k], '-', label=ZbinName[k], color=colors[k])
      plt.title(BoxName[i])

plt.legend()
plt.tight_layout()


fig = plt.figure(2, facecolor='w')
idx2 = np.nonzero(latCell<-60.0/180.0*pii)[0]
plt.plot(yCell[idx2], xCell[idx2], 'k.')
for i in range(nBoxes):
   idx = Boxes[i]
   plt.plot(yCell[idx], xCell[idx], '.', color=colors[i])
plt.show()


