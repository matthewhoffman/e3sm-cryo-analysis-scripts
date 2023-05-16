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
year_range = [150,200]

wed.tseries1(run_incr, var_incr,
             year_range=year_range, year_overlay=False, 
             placename = [placename],
             print_to_file=True, create_figure=False,
             #input_filename=filename,
             #apply_filter=True, cutoff=1/4,
             ztop_pyc=[True,False], zbottom_pyc=[False,True], 
             #diff_pyc=[True],
             #reference_run='ISMF', 
             savepath=savepath_anvil, overwrite=False)
             #zeval = [-100,-400],
             #zrange = [-100,-500],
             #lat=lat, lon=lon,
             #option = 'coord',placename = 'S4E',#'M31W',
             #zrange=[0,20],zab=True,

