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
from plot_config import M31W_T_obs, M31W_S_obs, M31W_rho_obs
run_incr = ['ISMF','ISMF-noEAIS','ISMF-3dGM','ISMF-noDIB']
placename = 'M31W'
year_range = [20,150]

#var_incr = ['T']
#wed.tseries1(run_incr,var_incr,year_range=year_range,
#             input_filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_M31W_z20.000000mab_t010-180',
#             placename = placename,zeval=20.,zab=True,
#             print_to_file=False,create_figure=True,
#             show_obs = True, obs = M31W_T_obs,
#             year_overlay=False,overwrite=False)
#
var_incr = ['S']
wed.tseries1(run_incr,var_incr,year_range=year_range,
             input_filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_M31W_z20.000000mab_t010-180',
             placename = placename,zeval=20.,zab=True,
             print_to_file=False,create_figure=True,
             show_obs = True, obs = M31W_S_obs,
             year_overlay=False,overwrite=False)

var_incr = ['rho']
wed.tseries1(run_incr,var_incr,year_range=year_range,
             input_filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_M31W_z20.000000mab_t010-180',
             placename = placename,zeval=20.,zab=True,
             print_to_file=False,create_figure=True,
             show_obs = True, obs = M31W_rho_obs,
             year_overlay=False,overwrite=False)

