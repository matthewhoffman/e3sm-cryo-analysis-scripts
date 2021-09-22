#!/usr/bin/env python
'''
Script to plot melt time series from different runs
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import csv
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

savepath=savepath_nersc
filename = 'ISMF_melt_t_Filchner_Ronne_log'
 
runs = OrderedDict()
runs['ISMF-noDIB'] = {'ts_dir':
#'/global/cscratch1/sd/dcomeau/e3sm_scratch/cori-knl/mpas-analysis-output/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/yrs120-150/timeseries/iceShelfFluxes',
'/global/homes/h/hoffman2/cryo_data/ISMF-noDIB',
'start':1,
'end':200
}

region=5 #FIS
Rregion=6 #RIS
FRregion=7#FRIS

#colors=[ cm.tab10(x) for x in np.linspace(0.0, 1.0, 4) ]
r=0
for run in runs:
    table_file = open(savepath + run + '_iceShelfMelt.txt','a+')
    wr = csv.writer(table_file,dialect='excel')
    wr.writerow(['decyear','melt_FIS','melt_RIS','melt_FRIS'])
    thisDict = runs[run]
    print(thisDict['ts_dir'])
    filelist = glob.glob(os.path.join(thisDict['ts_dir'], 'iceShelfFluxes*.nc'))

    Fmelt = np.zeros((thisDict['end']-thisDict['start']+1))
    Rmelt = np.zeros((thisDict['end']-thisDict['start']+1))
    FRmelt = np.zeros((thisDict['end']-thisDict['start']+1))
    year = np.zeros((thisDict['end']-thisDict['start']+1))
    i=0
    for infile in sorted(filelist): 
       #print(infile, i)
       #do some fancy stuff
       f=netCDF4.Dataset(infile, 'r')
       Fmelt[i] = f.variables['totalMeltFlux'][:, region].mean()
       year[i] = f.variables['Time'][0] / 365.0 + 1.0 # adding 1.0 offset because sims started then

       Rmelt[i] = f.variables['totalMeltFlux'][:, Rregion].mean()
       FRmelt[i] = f.variables['totalMeltFlux'][:, FRregion].mean()
       rowentries = [year[i], Fmelt[i], Rmelt[i], FRmelt[i]]
       wr.writerow(rowentries)
       i += 1
    r+=1
