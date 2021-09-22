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
run_incr = ['B-ISMF-bug','B-ISMF']
#year_range = [10,11]
year_range = [10,100]
savepath = '/home/ac.cbegeman/E3SM-analysis/weddell-fris/'
locname = 'trough_shelf'
filename = '_transect_u_'+locname+'_10-100'
#wed.fluxgate(locname,yrrange=year_range,morange = [1,12+1],run_incr=run_incr,
#             mode = 'barotropic-baroclinic', #mode = 'pos-neg',
#             plot_map = False, plot_transect = False,
#             overwrite=False,savepath = savepath)
#wed.tseries1(run_incr,
#             ['F_barotropic','F_baroclinic_pos'],
#             year_range=year_range,
#             placename=locname,
#             input_filename = run_incr[0] + filename,
#             savepath = savepath)
placename = 'M31W'
filename = '_S_M31W_z20.000000mab_t010-100'
wed.tseries1(run_incr,['T','S'],year_range=year_range,
             placename = placename,zeval=20.,zab=True,
             input_filename = [run_incr[0]+filename,run_incr[1]+filename],
             print_to_file=False,create_figure=True,
             savepath=savepath,year_overlay=False,overwrite=False)
#filename = '_rho_wed_pyc_Ryan_shelf_belowpyc_t010-100'
#placename = 'wed_pyc_Ryan_shelf'
#wed.tseries1(run_incr,['rho','rho'],year_range=year_range,
#             placename = placename,
#             input_filename = [run_incr[0]+filename,run_incr[1]+filename],
#             print_to_file=False,create_figure=True,
#             ztop_pyc = [True,False],
#             zbottom_pyc = [False,True],
#             savepath=savepath,year_overlay=False,overwrite=False)
