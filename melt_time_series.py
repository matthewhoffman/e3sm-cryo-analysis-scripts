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
from plot_config import *
from weddell_mod import extract_tseries

savepath=savepath_nersc
filename = 'ISMF_melt_t_Filchner_Ronne_log'

#run_list = ['ISMF-noDIB', 'ISMF', 'ISMF-3dGM']
run_list = ['ISMF-noDIB']
year_range = [0, 1]

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


runs['ISMF-3dGM'] = {'ts_dir':
#'/global/cscratch1/sd/xylar/analysis/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/',
'/global/homes/h/hoffman2/cryo_data/ISMF-MGM',
'start':21,
'end':100
}

Rignot_data = 'Ocean/Melt/Rignot_2013_melt_rates.csv'
# process obs
obsFileNameDict = {'Rignot et al. (2013)':
                   'Rignot_2013_melt_rates_20201117.csv',
                   'Rignot et al. (2013) SS':
                   'Rignot_2013_melt_rates_SS_20201117.csv'}

ice_shelves_to_plot = ['Filchner', 'Ronne']

# former is Rignot, latter is MPAS
ice_shelf_names = {'Filchner':'Filchner',
                   'Ronne':'Ronne',
                   'Brunt_Stancomb':'Brunt_Stancomb',
                   'Stancomb-Wills':None,
                   'Riiser-Larsen':None,
                   'Quar':None,
                   'Ekstrom':None
                  }
region=5 #FIS
Rregion=6 #RIS
#7=FRIS

extract_tseries(run_list, ['land_ice_fw_flux'], year_range,
                placename = 'EAcoast', land_ice_mask = True)

obsDict = {}  # dict for storing dict of obs data
# for ice_shelf_to_plot in ice_shelves_to_plot:
for obsName in obsFileNameDict:
    obsFileName = '{}/{}'.format(observationsDirectory,
                                 obsFileNameDict[obsName])
    obsDict[obsName] = {}
    obsFile = csv.reader(open(obsFileName, 'rU'))
    next(obsFile, None)  # skip the header line
    for line in obsFile:  # some later useful values commented out
        shelfName = line[0]
        if shelfName != ice_shelf_to_plot:
            continue

        # surveyArea = line[1]
        meltFlux = float(line[2])
        meltFluxUncertainty = float(line[3])
        meltRate = float(line[4])
        meltRateUncertainty = float(line[5])
        # actualArea = float( line[6] )  # actual area here is in sq km

        # build dict of obs. keyed to filename description
        # (which will be used for plotting)
        obsDict[obsName] = {
            'meltFlux': meltFlux,
            'meltFluxUncertainty': meltFluxUncertainty,
            'meltRate': meltRate,
            'meltRateUncertainty': meltRateUncertainty}
        break


fig1 = plt.figure(1, facecolor='w')
nrow=1
ncol=1
axMelt = fig1.add_subplot(nrow, ncol, 1)

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
    
    axMelt.plot(year, Fmelt, label=run, color=run_color[runname.index(run)],lineWidth=lw1)
    axMelt.plot(year, Rmelt, '--', color=run_color[runname.index(run)],lineWidth=lw1)
    r+=1


#plt.plot(ds.Time, ds.totalMeltFlux[:,0])
41.9
plt.xlabel('Year')
plt.ylabel('Ice shelf melt flux (Gt yr$^{-1}$)')
axMelt.plot([20,150], [41.9, 41.9], '-k', alpha=0.5, lineWidth=lw1)
axMelt.plot([20,150], [96.9, 96.9], '--k', alpha=0.5, lineWidth=lw1)
axMelt.set_yscale('log')
#plt.legend()

yl=axMelt.get_ylim()
#plt.plot([1, 1], [-100, 2000], 'k--', lineWidth=lw1)
#plt.plot([63, 63], [-100, 2000], 'k--', lineWidth=lw1)
#plt.plot([63+62, 63+62], [-100, 2000], 'k--', lineWidth=lw1)
for run in runs:
   plt.plot([run_tipping_year[runname.index(run)], run_tipping_year[runname.index(run)]], 
            [-100, 2000], ':', color=run_color[runname.index(run)],lineWidth=lw1)
axMelt.set_ylim(yl)
axMelt.set_xlim([20,150])

#plt.show()
plt.savefig(savepath + filename + '.png')#,dpi=set_dpi)
