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
from region_config import M31W_T_obs, M31W_S_obs
from data_config import runname

run_list = runname
placename = ['M31W']
year_range = [0, 200]
input_filename = ['CGM-UIB_T_S_M31W_z20.000000mab_t000-200',
                 'ISMF_ISMF-3dGM_T_S_M31W_z20mab_t001-200',
                 'ISMF_ISMF-3dGM_T_S_M31W_z20mab_t001-200']
zab = [True for i in run_list]
zeval = [20. for i in run_list]

# Compute tseries only for CGM-UIB, others already computed
#wed.tseries1(['CGM-UIB'], ['T','S'], year_range=year_range,
#             placename=placename, zeval=zeval, zab=zab,
#             print_to_file=True, create_figure=False,
#             overwrite=False)

wed.tseries1(run_list, ['T'], year_range=year_range,
             input_filename=input_filename,
             placename=placename, zeval=zeval, zab=zab,
             print_to_file=False, create_figure=True,
             show_obs='fill', obs=M31W_T_obs,
             year_overlay=False, overwrite=False)

wed.tseries1(run_list, ['S'], year_range=year_range,
             input_filename=input_filename,
             placename=placename, zeval=zeval, zab=zab,
             print_to_file=False, create_figure=True,
             show_legend=False, show_obs='fill', obs=M31W_S_obs,
             year_overlay=False, overwrite=False)
