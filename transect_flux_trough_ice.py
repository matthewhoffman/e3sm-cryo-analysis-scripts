#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm
import cmocean
from math import pi
import weddell_mod as wed
from pick_from_mesh import pick_transect 

#run_incr = ['ISMF','ISMF-noDIB','ISMF-noEAIS','ISMF-3dGM']
#run_incr = ['B-ISMF-bug']
run_incr = ['B-ISMF']
year_range = [10,11]
#year_range = [10,180]
savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/'
locname = 'trough_shelf'
#filename = '_transect_u_trough_ice_10-11_z.txt'
wed.fluxgate(locname,yrrange=year_range,morange = [1,12+1],run_incr=run_incr,
             mode = 'barotropic-baroclinic', #mode = 'pos-neg',
             plot_map = False, plot_transect = False,
             overwrite=False,savepath = savepath)
#for i in run_incr:
#    wed.hovmoller([i],year_range,option='transect',transect_id=locname,
#             input_filename = i + filename,
#             savepath = savepath,varlist = ['u_baroclinic'],limTrue = True)
