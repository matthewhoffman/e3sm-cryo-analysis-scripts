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

run_incr = ['ISMF','ISMF-noEAIS','ISMF-3dGM','ISMF-noDIB']
#placename = 'wed_pyc_Ryan_shallow'
#placename = 'gyre_interior'
placename = 'wed_pyc_west'
#var_incr = ['taux','tauy']
#var_incr = ['T','S','rho','u','v']
#var_incr = ['T','S']
#var_incr = ['unormal']
var_incr = ['rho']
rholim = [1027.2,1027.9]
lat = -70
lon = 340
#filename = 'ISMF_ISMF-noEAIS_rho_wed_pyc_Ryan_shallow_abovepyc_t070-101'
#filename = 'ISMF_zpyc_wed_pyc_west_10-180_zlim-800-10'
#filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_wed_pyc_west_abovepyc_t010-180'
filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_wed_pyc_west_belowpyc_t010-180'
year_range = [10,180]

#wed.tseries1(run_incr,['mean'],year_range=year_range,
#             placename = 'wed_pyc_Ryan', 
#             apply_filter = True, cutoff = 1/4,
#             print_to_file=True,create_figure=True,
#             input_filename = 'ISMF_zpyc_wed_pyc_Ryan_70-101_zlim-4500--3500')
#wed.tseries1(run_incr,['mean'],year_range=year_range,
#             placename = 'wed_pyc_Ryan', 
#             apply_filter = True, cutoff = 1/4,
#             print_to_file=True,create_figure=True,
#             input_filename = 'ISMF_zpyc_wed_pyc_Ryan_70-101_zlim-2000--500')
#wed.tseries1(run_incr,var_incr,year_range=year_range,
#             placename = 'gyre_interior',
#             print_to_file=True,create_figure=True,
#             apply_filter = True, cutoff = 1/4,
#             ztop_pyc = [True], zbottom_pyc = [False],
#             year_overlay=False,overwrite=True)
#wed.tseries1(run_incr,var_incr,year_range=year_range,
#             placename = placename,
#             print_to_file=True,create_figure=True,
#             varlim = rholim,
#             #apply_filter = True, cutoff = 1/4,
#             ztop_pyc = [False], zbottom_pyc = [True],
#             year_overlay=False,overwrite=True)
wed.tseries1(run_incr,var_incr,year_range=year_range,
             placename = placename,
             #reference_run='ISMF',
             #print_to_file=True,create_figure=False,
             print_to_file=False,create_figure=True,
             input_filename = filename,
             #apply_filter = True, cutoff = 1/4,
             #ztop_pyc = [True], zbottom_pyc = [False],
             ztop_pyc = [False], zbottom_pyc = [True],
             year_overlay=False,overwrite=True)
             #zeval = [-100,-400],
             #zrange = [-100,-500],
             #lat=lat, lon=lon,
             #option = 'coord',placename = 'S4E',#'M31W',
             #zrange=[0,20],zab=True,

