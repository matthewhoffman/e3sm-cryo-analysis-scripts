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
from plot_config import *

varlim = True
yr_incr = np.arange(191,201,1)
mo_incr = [6,7,8]
run_incr = ['ISMF','ISMF-3dGM']
#ops_incr = ['sigma1','']
var_incr = ['T']#,'S','sigma1','u','v']
pick_option = 'coord'
transect_name = [#'wtrough_crossshelf'
                 'trough_crossshelf'
                 #'etrough_crossshelf'
                ]
#levels = np.arange(1027.4,1028.0,0.1) # rho
levels = np.arange(32.0,32.7,0.05) # sigma

for tr in transect_name:
    for m in run_incr:
        wed.transect(pick_option, yr_incr, mo_incr, ['T'], 
                     transect_name = tr, tav=30,
                     overwrite = True, 
                     zlim = [5,-1e3],savepath=savepath_anvil,
                     run = m, plot_transect =False 
                    )
        #wed.transect(pick_option, yr_incr, mo_incr, ['S'], 
        #             transect_name = tr, tav=30,
        #             overwrite = True, 
        #             zlim = [5,-1e3],savepath=savepath_anvil,
        #             run = m, plot_transect =False 
        #            )
    #wed.transect(pick_option, [60], mo_incr, var_incr, 
    #             var_contour = 'sigma1',cntr_levels = levels,
    #             transect_name = tr, tav=3,
    #             overwrite = True, 
    #             zlim = [5,-1e3],savepath=savepath_anvil,
    #             run = 'ISMF-noDIB', plot_transect =False 
    #            )
