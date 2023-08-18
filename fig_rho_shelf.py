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
#from plot_config import gyre_interior_drho_obs, savepath_anvil

#placename = 'wed_pyc_Ryan_shelf'
placename = 'gyre_interior'
var_incr = ['rho','rho']
rholim = [1027.2, 1027.9]
filename = f'ISMF_ISMF-3dGM_rho_rho_{placename}_abovepyc_belowpyc_t001-200'#.txt'
filename_DIB=f'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_{placename}_abovepyc_t010-180'#.txt'
print(filename)
year_range = [25,200]

# Create timeseries
#wed.tseries1(['ISMF','ISMF-3dGM'], var_incr,
#             year_overlay=False, year_range=year_range,
#             placename=[placename],
#             print_to_file=True, create_figure=False,
#             input_filename=filename,
#             ztop_pyc=[True,False], zbottom_pyc=[False,True], 
#             savepath=savepath_anvil, overwrite=False)

# Plot timeseries
wed.tseries1(['CGM-UIB','CGM-DIB','VGM-DIB','CGM-UIB','CGM-DIB','VGM-DIB'], ['rho'],
             year_overlay=False, year_range=year_range,
             placename=[placename], varlim=False,
             apply_filter=True, cutoff=1/4,
             year_minmax=True, 
             linestyle=['--','--','--','-','-','-'],
             print_to_file=False, create_figure=True,
             input_filename=[filename_DIB,filename,filename,filename_DIB,filename,filename],
             ztop_pyc=[False, False, False, True, True, True], 
             zbottom_pyc=[True, True, True, False, False, False],
             #ztop_pyc=[False],
             #zbottom_pyc=[False], 
             diff_pyc=[False], #reference_run='ISMF', 
             #show_obs=True, obs=gyre_interior_drho_obs,
             #savepath=savepath_anvil,
             overwrite=False)
