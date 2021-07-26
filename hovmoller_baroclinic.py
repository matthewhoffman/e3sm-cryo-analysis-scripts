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

run_incr = ['ISMF','ISMF-noEAIS','ISMF-noDIB','ISMF-3dGM']
year_range = [10,180]
filename = ''

transect = 'trough_ice'
wed.hovmoller(i,startyr,endyr,option='transect',transect_id=transect,
             input_filename = filename,
             savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/',
             varlist = ['u_baroclinic'],limTrue = True)

