#!/usr/bin/env python
'''
'''

from netCDF4 import Dataset
from matplotlib import pyplot as plt
import numpy as np


CGMruns = [
             {'value':0.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/20210730.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.DIBbugfix.anvil/yrs191-200/timeseries/iceShelfFluxes/iceShelfFluxes_0001-0200.nc'},
#             {'value':16.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/CGM-16myr/yrs141-170/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0170.nc'},
             {'value':8.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/CGM-8myr/yrs141-150/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0150.nc'},
             {'value':4.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/CGM-4myr/yrs141-149/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0149.nc'},
#             {'value':2.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/CGM-2myr/yrs141-210/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0210.nc'},
#             {'value':0.5758, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/CGM-0.5758myr/yrs141-150/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0150.nc'}
]

VGMruns = [
             {'value':0.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/20210901.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.VGM.DIBbugfix.anvil//yrs191-200/timeseries/iceShelfFluxes/iceShelfFluxes_0001-0200.nc'},
             {'value':16.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/VGM-16myr/yrs141-170/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0170.nc'},
             {'value':8.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/VGM-8myr/yrs141-170/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0170.nc'},
             {'value':4.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/VGM-4myr/yrs141-180/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0180.nc'},
             {'value':2.0, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/VGM-2myr/yrs141-210/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0210.nc'},
             {'value':0.5758, 'path':'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/FRIS-branch-runs_Oct2021/VGM-0.5758myr/yrs141-150/timeseries/iceShelfFluxes/iceShelfFluxes_0141-0150.nc'}
]

# Calculate uniformIB high melt regime value
fnameUIB='/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/20191003.GMPAS-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/yrs001-101/timeseries/iceShelfFluxes/iceShelfFluxes_0001-0101.nc'
f = Dataset(fnameUIB, 'r')
time = f.variables['Time'][:] / 365.0
ind = time>87.0

region = 5 # Filchner
melt = f.variables['meltRates'][:, region]
UIBminFilchner = melt[ind].min()
UIBmaxFilchner = melt[ind].max()

region = 6 # Ronne
melt = f.variables['meltRates'][:, region]
UIBminRonne = melt[ind].min()
UIBmaxRonne = melt[ind].max()
f.close()




melts = np.array([16.0, 8.0, 4.0, 2.0, 0.5758])
CGMyrsToTP = np.zeros((len(melts),))
VGMyrsToTP = np.zeros((len(melts),))


fig = plt.figure(1, facecolor='w', figsize=(8,10))
colors = {0.0:'black', 16.0:'firebrick', 8.0:'tab:orange', 4.0:'gold', 2.0:'yellowgreen', 0.5758:'blue', 1.0:'limegreen'}


ax = fig.add_subplot(2, 1, 1)
region = 5 # Filchner
#region = 7 # FR
meltThreshold = 0.5758*1.5
plt.plot([0, 300], [UIBminFilchner, UIBminFilchner], ':', color='gray')
plt.plot([0, 300], [UIBmaxFilchner, UIBmaxFilchner], ':', color='gray')

i=0
for run in CGMruns:
   f = Dataset(run['path'], 'r')
   melt = f.variables['meltRates'][:, region]
   time = f.variables['Time'][:] / 365.0
   plt.plot(time, melt, color = colors[run['value']], label='{} m/yr'.format(run['value']))

   if run['value'] > 0.0:
      aryInd = np.where(melts==run['value'])[0][0]
      found = np.nonzero(melt>meltThreshold)[0]
      if len(found) > 0:
         ind = found[0]
         CGMyrsToTP[aryInd] = np.interp(meltThreshold, melt[ind-1:ind+1], time[ind-1:ind+1]) - 141.0
      else:
         CGMyrsToTP[aryInd] = np.nan
      print(CGMyrsToTP[aryInd])
   f.close()

   i+=1

i=0
for run in VGMruns:
   f = Dataset(run['path'], 'r')
   melt = f.variables['meltRates'][:, region]
   time = f.variables['Time'][:] / 365.0
   plt.plot(time, melt, '--', color = colors[run['value']], label='VGM {} m/yr'.format(run['value']))

   if run['value'] > 0.0:
      aryInd = np.where(melts==run['value'])[0][0]
      found = np.nonzero(melt>meltThreshold)[0]
      if len(found) > 0:
         ind = found[0]
         VGMyrsToTP[aryInd] = np.interp(meltThreshold, melt[ind-1:ind+1], time[ind-1:ind+1]) - 141.0
      else:
         VGMyrsToTP[aryInd] = np.nan
      print(VGMyrsToTP[aryInd])
   f.close()

   i+=1


plt.xlim([140, 210])
plt.xlabel('Year')
plt.ylabel('Filchner melt rate (m/yr)')
plt.legend()


ax = fig.add_subplot(2, 1, 2)
region = 6 # Ronne
plt.plot([0, 300], [UIBminRonne, UIBminRonne], ':', color='gray')
plt.plot([0, 300], [UIBmaxRonne, UIBmaxRonne], ':', color='gray')
i=0
for run in CGMruns:
   f = Dataset(run['path'], 'r')
   melt = f.variables['meltRates'][:, region]
   time = f.variables['Time'][:] / 365.0
   plt.plot(time, melt, color = colors[run['value']], label='{} m/yr'.format(run['value']))
   f.close()
   i+=1

i=0
for run in VGMruns:
   f = Dataset(run['path'], 'r')
   melt = f.variables['meltRates'][:, region]
   time = f.variables['Time'][:] / 365.0
   plt.plot(time, melt, '--', color = colors[run['value']])#, label='VGM {} m/yr'.format(run['value']))
   f.close()
   i+=1

plt.xlim([140, 210])
plt.xlabel('Year')
plt.ylabel('Ronne melt rate (m/yr)')


fig2 = plt.figure(2, facecolor='w')#, figsize=(8,10))
ax = fig2.add_subplot(1, 1, 1)
plt.plot(melts, CGMyrsToTP, 'ko', label='CGM')
for i in range(len(melts)):
   if np.isnan(CGMyrsToTP[i]):
      plt.plot(melts[i], 62.0, 'ro')
plt.plot(melts, VGMyrsToTP, 'kx', label='VGM')
for i in range(len(melts)):
   if np.isnan(VGMyrsToTP[i]):
      plt.plot(melts[i], 62.0, 'rx')
plt.legend()
plt.plot([-1.0, 17.0], [62.0, 62.0], 'r--')
plt.xlim([-0.2, 16.2])
plt.xlabel('Prescribed Eastern Weddell melt rate (m/yr)')
plt.ylabel('Years to tipping point')
plt.ylim([0, 65])

plt.show()
