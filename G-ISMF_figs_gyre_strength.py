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
from data_config import runname
from region_config import wedwang_c_btr_flux_obs

placename = ['wedwang_c']
var_incr = ['F_barotropic']#,'F_barotropic','F_barotropic']
filename = ['ISMF-noDIB_transect_u_wedwang_c_10-180',
            'ISMF_transect_u_wedwang_c_1-200',
            'ISMF_transect_u_wedwang_c_1-200'
            ]
year_range = [25, 200]

#wed.fluxgate(placename[0], run_incr=run_incr,
#             yrrange=year_range, morange = [1,12+1],
#             mode='barotropic-baroclinic', #mode = 'pos-neg',
#             plot_map = False, plot_transect=False,
#             overwrite=True, savepath=savepath_anvil)
wed.tseries1(runname, var_incr, year_range=year_range,
             placename=placename,
             #apply_filter = True, cutoff = 1/4,
             print_to_file=False, create_figure=True,
             input_filename=filename,
             varlim=True, flip_y=True, show_legend=True,
             show_obs=True, obs=wedwang_c_btr_flux_obs,
             #year_minmax=True, 
             overwrite=False)
             #ratio_barotropic = [True],
             #zeval = [-100,-400],
             #zrange = [-100,-500],
             #lat=lat, lon=lon,
             #option = 'coord',placename = 'S4E',#'M31W',
             #zrange=[0,20],zab=True,

