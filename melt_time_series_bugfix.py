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

savepath=savepath_anvil
filename = 'ISMF_melt_t_4panel_bugfix'
 
runs = OrderedDict()
runs['ISMF'] = {'ts_dir':
'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/20210730.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.DIBbugfix.anvil/yrs191-200/timeseries/iceShelfFluxes/',
#'/global/cscratch1/sd/dcomeau/e3sm_scratch/cori-knl/mpas-analysis-output/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/yrs120-150/timeseries/iceShelfFluxes',
#'/global/homes/h/hoffman2/cryo_data/ISMF',
'start':1,
'end':200
}

runs['ISMF-noDIB'] = {'ts_dir':
'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/20191003.GMPAS-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/yrs091-100/timeseries/iceShelfFluxes/',
#'/global/cscratch1/sd/dcomeau/e3sm_scratch/cori-knl/mpas-analysis-output/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/yrs120-150/timeseries/iceShelfFluxes',
#'/global/homes/h/hoffman2/cryo_data/ISMF-noDIB',
'start':1,
'end':100
}


runs['ISMF-3dGM'] = {'ts_dir':
'/lcrc/group/acme/ac.mhoffman/scratch/anvil/mpas_analysis_output/20210901.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.VGM.DIBbugfix.anvil/yrs191-200/timeseries/iceShelfFluxes/',
#'/global/cscratch1/sd/xylar/analysis/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/',
#'/global/homes/h/hoffman2/cryo_data/ISMF-MGM',
'start':1,
'end':200
}

fig1 = plt.figure(1, facecolor='w', figsize=(4,8))
nrow=4
ncol=1

regionNames = ["Antarctica", 
              "Peninsula", 
              "West Antarctica", 
              "East Antarctica", 
              "Larsen_C", 
              "Filchner", 
              "Ronne", 
              "Filchner-Ronne", 
              "Brunt-Stancomb", 
              "Fimbul", 
              "Amery", 
              "Totten", 
              "Eastern_Ross", 
              "Western_Ross", 
              "Ross", 
              "Getz", 
              "Thwaites", 
              "Pine_Island", 
              "Abbot", 
              "George_VI"
             ]
regionIDs = ["Antarctica", 
             #"Peninsula", 
             #"West Antarctica", 
             #"East Antarctica", 
             #"Larsen_C", 
             #"Filchner", 
             #"Ronne", 
             "Filchner-Ronne", 
             "Brunt-Stancomb", 
             "Fimbul", 
             #"Amery", 
             #"Totten", 
             #"Eastern_Ross", 
             #"Western_Ross", 
             #"Ross", 
             #"Getz", 
             #"Thwaites", 
             #"Pine_Island", 
             #"Abbot", 
             #"George_VI"
            ]
#colors=[ cm.tab10(x) for x in np.linspace(0.0, 1.0, 4) ]
for n,name in enumerate(regionIDs):
   axMelt = fig1.add_subplot(nrow, ncol, n+1)
   region = regionNames.index(name)
   print(name, region)
   for r,run in enumerate(runs):
       thisDict = runs[run]
       print(thisDict['ts_dir'])
       filelist = glob.glob(os.path.join(thisDict['ts_dir'], 'iceShelfFluxes*.nc'))
       data = np.zeros((thisDict['end']-thisDict['start']+2))
       year = np.zeros((thisDict['end']-thisDict['start']+2))

       i=0
       for infile in sorted(filelist): 
          f=netCDF4.Dataset(infile, 'r')
          year[i] = f.variables['Time'][0] / 365.0 + 1.0 # adding 1.0 offset because sims started then
          #print(infile,year[i])
          data[i] = f.variables['totalMeltFlux'][:, region].mean()
          i += 1
       
       axMelt.plot(year, data, label=runtitle[runname.index(run)], color=run_color[runname.index(run)],lineWidth=lw1)
       if n == 3:
          axMelt.set_xlabel('Year', fontsize=10)
       else:
          axMelt.set_xticklabels('')
       if n == 0:
          axMelt.legend(frameon=False, fontsize=10, loc='upper left')
       axMelt.set_ylabel('Melt flux (Gt yr$^{-1}$)', fontsize=10)
       axMelt.set_title(name, fontsize=10)
       #axMelt.set_yscale('log')

#plt.plot(ds.Time, ds.totalMeltFlux[:,0])
#axMelt.plot([20,150], [41.9, 41.9], '-k', alpha=0.5, lineWidth=lw1)
#axMelt.plot([20,150], [96.9, 96.9], '--k', alpha=0.5, lineWidth=lw1)
#fig1.tight_layout(pad=3.0)
#plt.legend()

yl=axMelt.get_ylim()
#plt.plot([1, 1], [-100, 2000], 'k--', lineWidth=lw1)
#plt.plot([63, 63], [-100, 2000], 'k--', lineWidth=lw1)
#plt.plot([63+62, 63+62], [-100, 2000], 'k--', lineWidth=lw1)
#for run in runs:
#   plt.plot([run_tipping_year[runname.index(run)], run_tipping_year[runname.index(run)]], 
#            [-100, 2000], ':', color=run_color[runname.index(run)],lineWidth=lw1)
#axMelt.set_ylim(yl)
#axMelt.set_xlim([20,150])

#plt.show()
plt.savefig(savepath + filename + '.png')#,dpi=set_dpi)
