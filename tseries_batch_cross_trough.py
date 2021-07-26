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
#placename = 'S4E'
#placename = 'M31W'
placename = 'N31W'
#var_incr = ['T','S','rho']#,'u','v']
#var_incr = ['T']#,'u','v']
var_incr = ['rho']#,'u','v']
#filename = 'ISMF_ISMF-noEAIS_rho_wed_pyc_Ryan_shallow_abovepyc_t070-101'
year_range = [10,180]
#year_range = [10,11]

wed.tseries1(run_incr,var_incr,year_range=year_range,
             #input_filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_S4E_z20.000000mab_t010-180',
             input_filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_M31W_z20.000000mab_t010-180',
             input_filename2 = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_N31W_z500.000000m_t010-180',
             #input_filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_rho_M31W_W_z20.000000mab_t010-180',
             var2 = 'rho_z500',
             placename = placename,zeval=20.,zab=True,
             #placename = placename,zeval=500.,zab=False,
             print_to_file=False,create_figure=True,
             #print_to_file=True,create_figure=False,
             year_overlay=False,overwrite=False)

