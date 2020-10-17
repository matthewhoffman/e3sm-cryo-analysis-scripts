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

run_incr = ['ISMF','ISMF-noEAIS']#,'ISMF-3dGM']
#lat = -1*(76 + 2.75/60)
#lon = 360 - (30 + 59.65/60)
#latS = 74.62
#lonW = 30.58
#latS = 76 + 2.75/60
#lonW = 30 + 59.65/60
lat = -73
lon = 330

for i in run_incr:

    wed.hovmoller(i,70,71,
                  option='coord',coord=[lat,lon],zlim=[0,-500],
                  varlist = ['T','S','rho'],limTrue=True,plot_pycnocline=True,
                  savepath = '/global/homes/c/cbegeman/weddell_output/')

#    startyr = 70
#    endyr = 80
#    transect = 'trough_shelf'
#    wed.hovmoller(i,startyr,endyr,option='transect',transect_id=transect,
#                 input_filename = i + '_transect_u_'+transect+'_' + str(startyr) + '-100.txt',
#                 savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/',
#                 varlist = ['u_barotropic','u_baroclinic'],limTrue = True)
#    transect = 'trough_ice'
#    wed.hovmoller(i,startyr,endyr,option='transect',transect_id=transect,
#                 input_filename = i + '_transect_u_'+transect+'_' + str(startyr) + '-100.txt',
#                 savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/',
#                 varlist = ['u_barotropic','u_baroclinic'],limTrue = True)
#    startyr = 90
#    endyr = 100
#    transect = 'trough_shelf'
#    wed.hovmoller(i,startyr,endyr,option='transect',transect_id=transect,
#                 input_filename = i + '_transect_u_'+transect+'_' + '70-100.txt',
#                 savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/',
#                 varlist = ['u_barotropic','u_baroclinic'],limTrue = True)
#    transect = 'trough_ice'
#    wed.hovmoller(i,startyr,endyr,option='transect',transect_id=transect,
#                 input_filename = i + '_transect_u_'+transect+'_' + '70-100.txt',
#                 savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/',
#                 varlist = ['u_barotropic','u_baroclinic'],limTrue = True)


