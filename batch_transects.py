#!/usr/bin/env python

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
from pick_from_mesh import *

# Directly east of the trough
lat_trough_east = [-76,-73]
lon_trough_east = [360-30,360-30]
#lat_trough_near_east = [-76,-73]
#lon_trough_near_east = [360-28,360-29]
lat_trough_near_east = [-75.5,-73]
lon_trough_near_east = [360-26,360-29]
lat_trough_mid_east = [-75.5,-72.5]
lon_trough_mid_east = [360-23,360-27]
lat_trough_far_east = [-73,-71.5]
lon_trough_far_east = [360-16,360-20]
lat_plot = lat_trough_near_east
lon_plot = lon_trough_near_east

varlim = True
yr_incr = np.arange(70,101,1)
#mo_incr = [2,5,8,11]
#mo_incr = [1,7]
mo_incr = np.arange(1,13,1)
run_incr = ['ISMF']#,'ISMF-noEAIS']#,'ISMF-3dGM']
#var_incr = ['rho']
#var_incr = ['u']
#ops_incr = ['sigma1','']
var_incr = ['T']
#var_incr = ['T','S','rho']#,'u','v']
#pick_option = 'by_index'
pick_option = 'coord'
overwrite = False
out = '/global/homes/c/cbegeman/weddell_output/fluxgate/'
transect_name = [#'wtrough_crossshelf'
                 'trough_crossshelf'
#                 ,'etrough_crossshelf'
                ]
#transect_name = ['trough_shelf']
#levels = np.arange(1027.4,1028.0,0.1)
levels = np.arange(32.0,32.7,0.05)
#levels = np.arange(-0.04,0.0,0.01)
#cellidx,edgeidx,dist,angle = pick_transect(option = 'by_index',
#                                   transect_name = 'trough_ice',
#                                   overwrite = False) 
#cellidx,edgeidx,dist,angle = pick_transect(option = 'coord',
#                                   lat = lat_plot, lon = lon_plot,
#                                   overwrite = True) 

for tr in transect_name:
    wed.transect(pick_option, [50], mo_incr, 
                 var_incr, #ops = ops_incr, 
                 var_contour = 'sigma1',cntr_levels = levels,
                 transect_name = tr,#lat = lat_plot, lon = lon_plot,
                 #ISW_contour = True,#zpyc_contour=True,
                 #normal= True,
                 overwrite = overwrite, 
                 run = run_incr[0], plot_transect=False
                 ,zlim = [5,-1e3]
                )
    for m in run_incr:
        wed.transect(pick_option, yr_incr, mo_incr, var_incr, 
                     var_contour = 'sigma1',cntr_levels = levels,
                     transect_name = tr,#lat = lat_plot, lon = lon_plot,
                     #ISW_contour = True,#zpyc_contour=True,
                     #normal= True,
                     overwrite = overwrite, 
                     run = m, plot_transect =False 
                     ,zlim = [5,-1e3]
                    )
