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
#idx=67250-1
#idx = np.argmin( (latCell-(-1.323514))**2 + (lonCell-5.672896)**2)  #122901-1
#idx=198673-1
idx = np.argmin( (latCell-(-77.75/180.0*pii))**2 + (lonCell-((360.0-36.15)/180.0*pii ))**2)  # Darelius 2016 Msouth
#idx = np.argmin( (latCell-(-77.0083/180.0*pii))**2 + (lonCell-((360.0-34.05)/180.0*pii ))**2)  # Darelius 2016 Mnorth
#idx=210384-1 # filchner sill
print "idx=",idx

maxLevelCell=fmesh.variables['maxLevelCell'][idx]
nz = maxLevelCell

bottomDepth = fmesh.variables['bottomDepth'][idx]
layerThickness = fmesh.variables['layerThickness'][0, idx, :nz]
thicknessSum = layerThickness.sum()
thicknessCumSum = layerThickness.cumsum()
zSurf = bottomDepth - thicknessSum
zLayerBot = zSurf - thicknessCumSum
z = zLayerBot + 0.5 * layerThickness

path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
#years = np.arange(70,105,3)
years = np.arange(60,80,1)
colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, len(years))]

fig = plt.figure(1, facecolor='w')
nrow=1
ncol=3

axT = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('temperature (deg. C)')
plt.ylabel('depth (m)')
plt.grid()

axS = fig.add_subplot(nrow, ncol, 2, sharey=axT)
plt.xlabel('salinity (psu)')
plt.ylabel('depth (m)')
plt.grid()

axrho = fig.add_subplot(nrow, ncol, 3, sharey=axT)
plt.xlabel('potential density (kg/m^3)')
plt.ylabel('depth (m)')
plt.grid()

for i in range(len(years)):
   yr = years[i]
   c = colors[i]
   print "yr=",yr
   f=netCDF4.Dataset(path+'/'+'mpaso.hist.am.timeSeriesStatsMonthly.{0:04d}-01-01.nc'.format(yr),'r')
   T=f.variables['timeMonthly_avg_activeTracers_temperature']
   S=f.variables['timeMonthly_avg_activeTracers_salinity']
   rho = f.variables['timeMonthly_avg_potentialDensity']

   axT.plot(T[0,idx,:nz], z, label="yr{0:04d}".format(yr), color=c)
   axT.plot(-0.0575*S[0,idx,:nz]+0.0901-7.61e-4*z, z, 'k:')
   axS.plot(S[0,idx,:nz], z, label="yr{0:04d}".format(yr), color=c)
   axrho.plot(rho[0,idx,:nz], z, label="yr{0:04d}".format(yr), color=c)

#axT.set_xlim([-10,20])
#axS.set_xlim([28,36])
#axrho.set_xlim([1000,1040])

#axT.legend()
axS.legend()

#axT.plot([-1.8, -1.8], [-4000.0, 0.0], 'k:')


plt.show()
