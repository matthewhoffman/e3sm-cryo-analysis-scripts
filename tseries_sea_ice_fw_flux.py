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

#run_incr = ['ISMF','ISMF-noEAIS','ISMF-noDIB']
run_incr = ['ISMF','ISMF-noEAIS','ISMF-3dGM','ISMF-noDIB']
placename = 'frisEAcoast'
var_incr = ['sea_ice_fw_flux']
year_range = [10,180]

wed.tseries1(run_incr,var_incr,year_range=year_range,
             placename = placename,
             input_filename = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_sea_ice_fw_flux_frisEAcoast_t010-180',
             print_to_file=False,create_figure=True,
             operation = 'area_sum',
             #reference_run = 'ISMF',
             #tav = 12,
             #apply_filter = True, cutoff = 1/4,
             year_overlay=True,overwrite=False)
