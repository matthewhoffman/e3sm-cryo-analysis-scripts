#!/usr/bin/env python
'''
Script to plot melt time series from different runs
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
from collections import OrderedDict
import xarray as xr
import glob


runs = OrderedDict()
runs['ISMF'] = {'ts_dir':
#'/global/cscratch1/sd/dcomeau/e3sm_scratch/cori-knl/mpas-analysis-output/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/yrs120-150/timeseries/iceShelfFluxes',
'/global/homes/h/hoffman2/cryo_data/ISMF',
'start':1,
'end':150
}

runs['ISMF-noEAIS'] = {'ts_dir':
#'/global/cscratch1/sd/dcomeau/e3sm_scratch/cori-knl/mpas-analysis-output/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/yrs120-150/timeseries/iceShelfFluxes',
'/global/homes/h/hoffman2/cryo_data/ISMF-noEAIS',
'start':50,
'end':150
}

runs['ISMF-noDIB'] = {'ts_dir':
#'/global/cscratch1/sd/dcomeau/e3sm_scratch/cori-knl/mpas-analysis-output/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/yrs120-150/timeseries/iceShelfFluxes',
'/global/homes/h/hoffman2/cryo_data/ISMF-noDIB',
'start':1,
'end':100
}


runs['ISMF-MGM'] = {'ts_dir':
#'/global/cscratch1/sd/xylar/analysis/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/',
'/global/homes/h/hoffman2/cryo_data/ISMF-MGM',
'start':21,
'end':100
}

fig1 = plt.figure(1, facecolor='w')
nrow=1
ncol=1
axMelt = fig1.add_subplot(nrow, ncol, 1)

region=5 #FIS
Rregion=6 #RIS
#7=FRIS

#colors=[ cm.tab10(x) for x in np.linspace(0.0, 1.0, 4) ]
r=0
for run in runs:
    thisDict = runs[run]
    print(thisDict['ts_dir'])
    filelist = glob.glob(os.path.join(thisDict['ts_dir'], 'iceShelfFluxes*.nc'))

    Fmelt = np.zeros((thisDict['end']-thisDict['start']+1))
    Rmelt = np.zeros((thisDict['end']-thisDict['start']+1))
    year = np.zeros((thisDict['end']-thisDict['start']+1))
    i=0
    for infile in sorted(filelist): 
       #print(infile, i)
       #do some fancy stuff
       f=netCDF4.Dataset(infile, 'r')
       Fmelt[i] = f.variables['totalMeltFlux'][:, region].mean()
       year[i] = f.variables['Time'][0] / 365.0 + 1.0 # adding 1.0 offset because sims started then

       Rmelt[i] = f.variables['totalMeltFlux'][:, Rregion].mean()
       i += 1
    
    axMelt.plot(year, Fmelt, label=run, color=cm.tab10(r))
    axMelt.plot(year, Rmelt, '--', color=cm.tab10(r))
    r+=1



#plt.plot(ds.Time, ds.totalMeltFlux[:,0])

plt.xlabel('Year')
plt.ylabel('Filchner Ice Shelf melt flux (Gt yr$^{-1}$)')
plt.legend()

yl=axMelt.get_ylim()
plt.plot([1, 1], [-100, 2000], 'k--')
plt.plot([63, 63], [-100, 2000], 'k--')
plt.plot([63+62, 63+62], [-100, 2000], 'k--')
axMelt.set_ylim(yl)

plt.show()
