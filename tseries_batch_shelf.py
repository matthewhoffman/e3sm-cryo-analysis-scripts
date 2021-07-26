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
placename = 'wed_pyc_Ryan_shelf'
#var_incr = ['taux','tauy']
#var_incr = ['T','S','rho','u','v']
#var_incr = ['T','S']
#var_incr = ['unormal']
var_incr = ['rho']
#filename = 'ISMF_ISMF-noEAIS_rho_wed_pyc_Ryan_shallow_abovepyc_t070-101'
year_range = [10,180]

wed.tseries1(run_incr,var_incr,year_range=year_range,
             placename = placename,
             print_to_file=False,create_figure=True,
             #input_filename = filename,
             #apply_filter = True, cutoff = 1/4,
             ztop_pyc = [True for i in run_incr],
             year_overlay=False,overwrite=False)
