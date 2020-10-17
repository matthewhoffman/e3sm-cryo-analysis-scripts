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
lat_plot = lat_trough_east
lon_plot = lon_trough_east

varlim = True
#yr_incr = np.arange(90,100,1)
yr_incr = [94]
mo_incr = [2,5,8,11]
#mo_incr = [1]
#mo_incr = np.arange(1,13,1)
#run_incr = ['ISMF','ISMF-noEAIS']#,'ISMF-3dGM']
run_incr = ['ISMF-noEAIS']#,'ISMF-3dGM']
#var_incr = ['rho']
var_incr = ['T','S','rho']#,'u','v']
out = '/global/homes/c/cbegeman/weddell_output/fluxgate/'

#cellidx,edgeidx,dist,angle = pick_transect(option = 'by_index',
#                                   transect_name = 'trough_ice',
#                                   overwrite = False) 
#cellidx,edgeidx,dist,angle = pick_transect(option = 'coord',
#                                   lat = lat_plot, lon = lon_plot,
#                                   overwrite = True) 

for m in run_incr:
    wed.transect('coord', yr_incr, mo_incr, var_incr, 
                 lat = lat_plot, lon = lon_plot,
                 varlim = varlim, overwrite = False, 
                 run = m, plot_transect =False 
                 ,zlim = [5,-1e3]
                )
    wed.transect('coord', yr_incr, mo_incr, ['u'], 
                 lat = lat_plot, lon = lon_plot,
                 varlim = varlim,overwrite=False, 
                 run = m,plot_transect=False
                 ,zlim = [5,-1e3],normal=True
                )
