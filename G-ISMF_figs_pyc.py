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
from math import pi
import weddell_mod as wed
from plot_config import savepath_anvil

run_incr = ['ISMF','ISMF-3dGM']
#placename = 'wed_pyc_Ryan_shallow'
placename = 'wed_pyc_Ryan_shelf'
var_incr = ['rho','rho']
rholim = [1027.2,1027.9]
#filename = '_'.join(run_incr)+'_rho_'+placename+'_abovepyc_t001-002'
filename='ISMF_ISMF-3dGM_rho_rho_wed_pyc_Ryan_shelf_abovepyc_belowpyc_t001-200.txt'
filename_noDIB='ISMF-noDIB_data/ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_wed_pyc_Ryan_shelf_abovepyc_t010-180.txt'
year_range = [10,200]

# Create timeseries
#wed.tseries1(run_incr, var_incr,
#             year_range=year_range, year_overlay=False, 
#             placename = [placename],
#             print_to_file=True, create_figure=False,
#             ztop_pyc=[True,False], zbottom_pyc=[False,True], 
#             savepath=savepath_anvil, overwrite=False)

# Plot timeseries
wed.tseries1(['ISMF-noDIB','ISMF-3dGM','ISMF'], ['rho'],
             year_range=year_range, 
             placename = [placename], varlim=True,
             print_to_file=False, create_figure=True,
             apply_filter = True, cutoff = 1/4,
             input_filename=[filename_noDIB,filename,filename],
             plot_legend=False,
             #ztop_pyc=[True], zbottom_pyc=[False], 
             year_overlay=False, year_minmax=True,
             ztop_pyc=[False], zbottom_pyc=[True], 
             savepath=savepath_anvil, overwrite=False)
