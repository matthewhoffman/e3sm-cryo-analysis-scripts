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

#fmesh=netCDF4.Dataset('/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
fmesh=netCDF4.Dataset('/lcrc/group/e3sm/data/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
latCell = fmesh.variables['latCell'][:]
lonCell = fmesh.variables['lonCell'][:]
xCell = fmesh.variables['xCell'][:]
yCell = fmesh.variables['yCell'][:]
depths = fmesh.variables['refBottomDepth'][:]
#idx=67250-1
#idx = np.argmin( (latCell-(-1.323514))**2 + (lonCell-5.672896)**2)  #122901-1
#idx=198673-1
#idx = np.argmin( (latCell-(-77.75/180.0*pii))**2 + (lonCell-((360.0-36.15)/180.0*pii ))**2)  # Darelius 2016 Msouth
#idx = np.argmin( (latCell-(-77.0083/180.0*pii))**2 + (lonCell-((360.0-34.05)/180.0*pii ))**2)  # Darelius 2016 Mnorth
#idx = np.argmin( (latCell-(-74.67/180.0*pii))**2 + (lonCell-((360.0-30.5)/180.0*pii ))**2)  # S4E Arthun et al 2012
idx=210384-1 # filchner sill
print("idx=",idx)

maxLevelCell=fmesh.variables['maxLevelCell'][idx]
nz = maxLevelCell

bottomDepth = fmesh.variables['bottomDepth'][idx]
layerThickness = fmesh.variables['layerThickness'][0, idx, :nz]
thicknessSum = layerThickness.sum()
thicknessCumSum = layerThickness.cumsum()
zSurf = bottomDepth - thicknessSum
zLayerBot = zSurf - thicknessCumSum
z = zLayerBot + 0.5 * layerThickness
print(z)



# plot location
fig = plt.figure(2, facecolor='w')
idx2 = np.nonzero(latCell<-60.0/180.0*pii)[0]
plt.plot(yCell[idx2], xCell[idx2], 'k.')
plt.plot(yCell[idx], xCell[idx], 'r.')

# OLD PATHS
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190923.GMPAS-IAF.T62_oEC60to30v3wLI.cori-knl/run'
#path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist'
#path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist'
##path='/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'

##path='/global/cscratch1/sd/kehoch/acme_scratch/cori-knl/20190304.GMPAS-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
##path='/global/cscratch1/sd/kehoch/acme_scratch/cori-knl/20190304.GMPAS-IAF.T62_oEC60to30v3wLI.cori-knl/run'

##path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'

##path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.edison/ocn/hist'
##path='/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190306.A_WCYCL1850-DIB_CMIP6.ne30_oECv3wLI.edison/ocn/hist'

#path='/global/cscratch1/sd/sprice/acme_scratch/cori-knl/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/run'

#path='/global/cscratch1/sd/sprice/acme_scratch/cori-knl/20191018.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.cori-knl.testNewGM/run'

# NEW PATHS in Project dir
#path='/project/projectdirs/m3412/simulations/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist/'
#path='/project/projectdirs/m3412/simulations/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist/'
path='/lcrc/group/acme/ac.mhoffman/acme_scratch/anvil/20210730.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.DIBbugfix.anvil/run'
path='/lcrc/group/acme/ac.mhoffman/acme_scratch/anvil/20210818.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.DIBbugfix-17myrDML.anvil/run'
path='/lcrc/group/acme/ac.dcomeau/acme_scratch/anvil/20210614.A_WCYCL1850-DIB-ISMF_CMIP6.ne30_oECv3wLI.DIBbugfix.anvil/run'
path='/lcrc/group/acme/ac.mhoffman/acme_scratch/anvil/20210901.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.VGM.DIBbugfix.anvil/run'

years = np.arange(50,129,1)
years = np.arange(1,139,1)
years = np.arange(1,10,1)
#years = np.arange(65,75,1)
#years = np.arange(1,100,1)
years = np.arange(1,139,1)
#years = np.arange(201,230,1)
months = np.arange(1,13,1)
nt = len(years)*len(months)
times = np.zeros((nt,))


Tdata = np.zeros((nz, nt))
Sdata = np.zeros((nz, nt))
rhodata = np.zeros((nz, nt))
vdata = np.zeros((nz, nt))
# FWF
SIdata = np.zeros((nt,))
PEdata = np.zeros((nt,))
AISdata = np.zeros((nt,))
Restoredata = np.zeros((nt,))
SSSData = np.zeros((nt,))
SSHData = np.zeros((nt,))
S100Data = np.zeros((nt,))
wsData = np.zeros((nt,))
wsDirData = np.zeros((nt,))
velData = np.zeros((nt,))
velDirData = np.zeros((nt,))
uwsData = np.zeros((nt,))
vwsData = np.zeros((nt,))
uData = np.zeros((nt,))
vData = np.zeros((nt,))

k100 = np.nonzero(z<=-80.0)[0][0]

t=0
for yr in years:
 for mo in months:
   print("yr=",yr, "mo=", mo)
   times[t] = yr+(mo-1.0)/12.0

   f=netCDF4.Dataset('{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo), 'r')
   T=f.variables['timeMonthly_avg_activeTracers_temperature']
   S=f.variables['timeMonthly_avg_activeTracers_salinity']
   rho = f.variables['timeMonthly_avg_potentialDensity']
   v  = f.variables['timeMonthly_avg_velocityMeridional']
   
   # store data
   Tdata[:,t] = T[0,idx,:nz]
   Sdata[:,t] = S[0,idx,:nz]
   rhodata[:,t] = rho[0,idx,:nz]
   vdata[:,t] = v[0,idx,:nz]

   # FWF data
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
   SSS=f.variables['timeMonthly_avg_activeTracers_salinity'][0,idx,0]
   dens=f.variables['timeMonthly_avg_potentialDensity'][0,idx,0]
   if 'timeMonthly_avg_salinitySurfaceRestoringTendency' in f.variables:
      restTend=f.variables['timeMonthly_avg_salinitySurfaceRestoringTendency'][0,idx]
   else:
       restTend = si*0.0
   SSH=f.variables['timeMonthly_avg_pressureAdjustedSSH'][0,idx]
   vws=f.variables['timeMonthly_avg_windStressZonal'][0,idx]
   uws=f.variables['timeMonthly_avg_windStressMeridional'][0,idx]
   u=f.variables['timeMonthly_avg_velocityMeridional'][0,idx,0]
   v=f.variables['timeMonthly_avg_velocityZonal'][0,idx,0]

   # store data
   SIdata[t] = si
   PEdata[t] = rain+snow+evap
   AISdata[t] = river + iceberg + iceshelf
   Restoredata[t] = -1.0*restTend/SSS*dens
   SSSData[t] = SSS
   S100Data[t] = S[0,idx,k100]
   SSHData[t] = SSH
   vwsData[t] = vws
   uwsData[t] = uws
   uData[t] = u
   vData[t] = v

   f.close()
   t += 1
 

fig = plt.figure(1, facecolor='w')
nrow=6
ncol=1

maxDepth=-200.0

axT = fig.add_subplot(nrow, ncol, 2)
plt.xlabel('year')
plt.ylabel('temperature\n(deg. C)')
plt.pcolor(times, z[:nz], Tdata, vmin=-2.1, vmax=1.6, cmap='nipy_spectral')#, vmin=-2.0, vmax=0.0)
#plt.pcolor(times, z[:nz], Tdata, vmin=-2.25, vmax=-1.5, cmap='RdBu_r')# Darelius 2016 approx.
#plt.pcolor(times, z[:nz], Tdata, vmin=-2.2, vmax=-0.5, cmap=plt.cm.get_cmap('inferno', 17))# Daae GRL figure
plt.pcolor(times, z[:nz], Tdata, vmin=-2.25, vmax=-1.40, cmap=plt.cm.get_cmap('RdYlBu_r', 17))# Darelius 2016  figure 2
#plt.pcolor(times, z[:nz], Tdata, vmin=-2., vmax=-1.0, cmap=plt.cm.get_cmap('RdYlBu_r', 17))# Arthun 2012  figure 3
#axT.set_ylim((maxDepth, 0))
plt.colorbar()
#plt.contour(times, z[:nz], Tdata, [-1.8])
plt.title(path + ", idx={}".format(idx))


axS = fig.add_subplot(nrow, ncol, 3, sharex=axT)
plt.xlabel('year')
plt.ylabel('salinity\n(psu)')
plt.pcolor(times, z[:nz], Sdata, vmin=33.5,vmax=34.6, cmap='nipy_spectral')#, vmin=31.0, vmax=34.7)
#axS.set_ylim((maxDepth, 0))
plt.colorbar()

axrho = fig.add_subplot(nrow, ncol, 4, sharex=axT)
plt.xlabel('year')
plt.ylabel('pot. dens.\n(kg/m3)')
plt.pcolor(times, z[:nz], rhodata, vmin=1027.2,vmax=1027.9, cmap='nipy_spectral')#, vmin=31.0, vmax=34.7)
#axS.set_ylim((maxDepth, 0))
plt.colorbar()

axv = fig.add_subplot(nrow, ncol, 5, sharex=axT)
plt.xlabel('year')
plt.ylabel('n velo\n(m/s)')
plt.pcolor(times, z[:nz], vdata, cmap='seismic', vmin=-0.05, vmax=0.05)
#axS.set_ylim((maxDepth, 0))
plt.colorbar()


waterMassMask = np.zeros((nz, nt))
waterMassMask[np.where((Sdata>34.0) * (Tdata>0.0))] = 1 # CDW
waterMassMask[np.where((Sdata>34.0) * (Tdata<0.0) * (Tdata>-1.5))] = 2 # mCDW
waterMassMask[np.where( (Tdata>(-0.0575*Sdata+0.0901)) * (Sdata>34.5) * (Tdata<-1.5))] = 3 # HSSW
waterMassMask[np.where( (Tdata>(-0.0575*Sdata+0.0901)) * (Sdata<34.5) * (Sdata>34.0) * (Tdata<-1.5))] = 4 # LSSW
waterMassMask[np.where( (Tdata>(-0.0575*Sdata+0.0901)) * (Sdata<34.0))] = 5 # AASW
waterMassMask[np.where(Tdata<(-0.0575*Sdata+0.0901))] = 6 # ISW

#fig = plt.figure(5, facecolor='w')
axMass = fig.add_subplot(nrow,ncol,6, sharex=axT)
plt.pcolor(times, z[:nz], waterMassMask, vmin=0.5, vmax=6.5, cmap=cm.get_cmap('tab10',6) )
cbar = plt.colorbar(ticks=[1,2,3,4,5,6])
cbar.set_ticklabels(['CDW','mCDW','HSSW','LSSW','AASW','ISW']) 
plt.xlabel('year')
plt.ylabel('depth (m)')

axFWF = fig.add_subplot(nrow, ncol, 1, sharex=axT)
plt.xlabel('year')
plt.ylabel('FWF (m/yr)')
factor = 1.0/1000.0*3.14e7
plt.plot(times, AISdata*factor, label='iceberg+shelf')
plt.plot(times, SIdata*factor, label='sea ice')
plt.plot(times, PEdata*factor, label='P-E')
plt.plot(times, Restoredata*factor, label='restoring')
plt.plot(times, (AISdata+SIdata+PEdata+Restoredata)*factor, 'k', label='Total')
#plt.plot(times, SSSData-34, label='sal-34')
plt.legend()
plt.colorbar()




fig = plt.figure(5, facecolor='w')
nrow=7
ncol=1

FWFdata=AISdata+SIdata+PEdata+Restoredata

doFilter = 1
if doFilter:
   filterLen=12
   uwsData = np.convolve(np.ones(filterLen,'d')/filterLen, uwsData, mode='same')    
   vwsData = np.convolve(np.ones(filterLen,'d')/filterLen, vwsData, mode='same')    
   uData = np.convolve(np.ones(filterLen,'d')/filterLen, uData, mode='same')    
   vData = np.convolve(np.ones(filterLen,'d')/filterLen, vData, mode='same')    
   #SSSData = np.convolve(np.ones(filterLen,'d')/filterLen, SSSData, mode='same')    
   #S100Data = np.convolve(np.ones(filterLen,'d')/filterLen, S100Data, mode='same')    
   FWFdata = np.convolve(np.ones(filterLen,'d')/filterLen, FWFdata, mode='same')    

wsData = (uwsData**2+vwsData**2)**0.5
wsDirData = np.arctan2(vwsData,uwsData)*180.0/3.14159
velData = (uData**2+vData**2)**0.5
velDirData = np.arctan2(vData,uData)*180.0/3.14159


axSal = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('year')
plt.ylabel('salinity')
plt.plot(times, SSSData, label='SSS')
plt.plot(times, S100Data, label='S80')
plt.legend()

axFWFTS = fig.add_subplot(nrow, ncol, 2, sharex=axSal)
plt.xlabel('year')
plt.ylabel('FWF')
plt.plot(times, (FWFdata)*factor, 'k', label='Total')

axSSH = fig.add_subplot(nrow, ncol, 3, sharex=axSal)
plt.xlabel('year')
plt.ylabel('SSH')
plt.plot(times, SSHData)

axWS = fig.add_subplot(nrow, ncol, 4, sharex=axSal)
plt.xlabel('year')
plt.ylabel('wind stress mag (N/m2)')
plt.plot(times, wsData)

axWS2 = fig.add_subplot(nrow, ncol, 5, sharex=axSal)
plt.xlabel('year')
plt.ylabel('wind stress dir ')
plt.plot(times, wsDirData)

axWS2 = fig.add_subplot(nrow, ncol, 6, sharex=axSal)
plt.xlabel('year')
plt.ylabel('surf vel mag (m/s) ')
plt.plot(times, velData)

axWS2 = fig.add_subplot(nrow, ncol, 7, sharex=axSal)
plt.xlabel('year')
plt.ylabel('surf vel dir (m/s) ')
plt.plot(times, velDirData)






plt.show()
