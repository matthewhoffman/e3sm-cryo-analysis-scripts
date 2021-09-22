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
from math import pi,nan
import weddell_mod as wed
from pick_from_mesh import pick_transect 

run_incr = ['ISMF','ISMF-noDIB','ISMF-3dGM']
year_range = [20,21]
#year_range = [10,100]
savepath = '/home/ac.cbegeman/E3SM-analysis/weddell-fris-bugfix/'
placename = 'M31W'

#wed.tseries1(run_incr,['land_ice_fw_flux','T','S'],year_range=year_range,
#             placename = ['',placename,placename],zeval=[-9999,20.,20.],zab=[True,True,True],
#             print_to_file=True,create_figure=True,
#             savepath=savepath,year_overlay=False,overwrite=False)
wed.tseries1(run_incr,['T','S'],year_range=year_range,
             placename = [placename,placename],zeval=[20.,20.],zab=[True,True],
             varlim = True, print_to_file=True,create_figure=True,
             savepath=savepath,year_overlay=False,overwrite=False)
