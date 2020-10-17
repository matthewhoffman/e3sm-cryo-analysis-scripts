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
#from transects import transect

latS = [82,73]
lonW = [70,0]
lonW_2 = [50,30]
lonW_3 = [48,30]
#lat_incr = np.arange(74,83,2)#1)
#lon_incr = np.arange(30,40,2)#1)
varlim = True
yr_incr = np.arange(80,81,1)
mo_incr = np.arange(1,2,1)
run_incr = ['ISMF','ISMF-noEAIS']
var_incr = ['T','S','rho','u','v']
out = '/global/homes/c/cbegeman/weddell_output/new_transects/'

for m in run_incr:
    for j in yr_incr:
        for k in mo_incr:
            wed.transect([76,76],lonW,var_incr,j,k,varlim,run=m,savepath=out)
            wed.transect([77,77],lonW,var_incr,j,k,varlim,run=m,savepath=out)
            wed.transect([78,78],lonW_3,var_incr,j,k,varlim,run=m,savepath=out)
            wed.transect([79,79],lonW_2,var_incr,j,k,varlim,run=m,savepath=out)
            wed.transect(latS,[38,38],var_incr,j,k,varlim,run=m,savepath=out)
            wed.transect(latS,[32,32],var_incr,j,k,varlim,run=m,savepath=out)
