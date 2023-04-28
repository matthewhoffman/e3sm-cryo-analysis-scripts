#!/usr/bin/env python
'''
Script to plot melt time series from different runs
'''

import os
import pandas
from data_config import *
from extract_data import extract_tseries
from plot_config import *

run_list = runname
year_range = [10, 100]

#extract_tseries([run_list], ['land_ice_fw_flux'], year_range,
#                placename = 'fris', land_ice_mask = False)
var = 'land_ice_fw_flux'
filename = f'{savepath}/bogus'
if not os.path.exists(filename):
    extract_tseries(run_list, [var], year_range,
                    placename='EAcoast', land_ice_mask=True,
                    operation='area_mean')

#i=0
#df = pandas.read_csv(filename)
#times = df['decyear'][:]
#header = f'{runtitle[i]}_{var}'
#melt_flux = df[header][:]
#
## the output is now in units of kg s^-1 m^-2
#melt_rate = melt_flux * s_to_yr / rho_ice
#
#print(melt_rate)
