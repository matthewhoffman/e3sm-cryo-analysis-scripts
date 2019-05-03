#!/usr/bin/env python

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm

fmesh=netCDF4.Dataset('oEC60to30v3wLI60lev.171031.nc')
#fmesh=netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]

path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'
years = np.arange(50,52,1)
months = np.arange(1,13,1)
nt = len(years)*len(months)
times = np.zeros((nt))

fw_flux = np.zeros((nt))

t=0
for yr in years:
 for mo in months:
   print "yr=",yr, "mo=", mo
   times[t] = yr+(mo-1.0)/12.0

   f=netCDF4.Dataset(path+'/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   temp = f.variables['timeMonthly_avg_seaIceFreshWaterFlux']
   
   # get data
   fw_flux[t] = np.sum(np.multiply(temp[0,:],fmesh.variables['areaCell'][:]))/np.sum(fmesh.variables['areaCell'][:])
   
   f.close()
   t += 1

fig = plt.figure()#1, facecolor='w')
plt.plot(times,fw_flux,'k')
plt.ylabel('Sea Ice Freshwater flux (kg m^(-2) s^(-1))')
plt.xlabel('Year')

plt.savefig('SeaIceFreshwaterFlux_globalavg_monthlyavg_tseries.png')
