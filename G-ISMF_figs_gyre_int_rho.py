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
from plot_config import gyre_interior_drho_obs, savepath_anvil

run_incr = ['ISMF','ISMF-3dGM']
placename = 'gyre_interior'
var_incr = ['rho','rho']
rholim = [1027.2,1027.9]
filename = 'ISMF_ISMF-3dGM_rho_rho_gyre_interior_abovepyc_belowpyc_t001-200.txt'
filename_DIB='ISMF-noDIB_data/ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_gyre_interior_abovepyc_t010-180.txt'
print(filename)
year_range = [10,200]

# Create timeseries
#wed.tseries1(['ISMF','ISMF-3dGM'], var_incr,
#             year_overlay=False, year_range=year_range,
#             placename=[placename],
#             print_to_file=True, create_figure=False,
#             input_filename=filename,
#             ztop_pyc=[True,False], zbottom_pyc=[False,True], 
#             savepath=savepath_anvil, overwrite=False)

# Plot timeseries
wed.tseries1(['ISMF-3dGM','ISMF','ISMF-noDIB'], ['rho'],
             year_overlay=False, year_range=year_range,
             placename=[placename], varlim=False,
             print_to_file=False, create_figure=True,
             input_filename=[filename,filename,filename_DIB],
             ztop_pyc=[True], zbottom_pyc=[False], 
             diff_pyc=[True], #reference_run='ISMF', 
             show_obs=True, obs=gyre_interior_drho_obs,
             savepath=savepath_anvil, overwrite=False)
