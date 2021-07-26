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
from plot_config import wedwang_c_btr_flux_obs

run_incr = ['ISMF','ISMF-noEAIS','ISMF-3dGM','ISMF-noDIB']
placename = 'wedwang_c'
#var_incr = ['F_baroclinic_pos']
var_incr = ['F_barotropic']
savepath_flux = '/global/homes/c/cbegeman/weddell_output/fluxgate/'
filename = 'ISMF_transect_u_wedwang_c_10-180'
year_range = [10,180]

wed.tseries1(run_incr,var_incr,year_range=year_range,
             placename = placename,
             print_to_file=False,create_figure=True,
             input_filename = filename,
             savepath = savepath_flux,
             show_obs = True, obs = wedwang_c_btr_flux_obs,
             year_overlay=False,overwrite=True)
             #ratio_barotropic = [True],
             #apply_filter = True, cutoff = 1/4,
             #zeval = [-100,-400],
             #zrange = [-100,-500],
             #lat=lat, lon=lon,
             #option = 'coord',placename = 'S4E',#'M31W',
             #zrange=[0,20],zab=True,

