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

lat_incr = np.arange(74,83,2)#1)
lon_incr = np.arange(30,40,2)#1)
yr = 70
mo = 1
varlim = True
yr_incr = np.arange(94,111,2)
mo_incr = np.arange(1,2,3)
run_incr = ['ISMF-noEAIS']
var_incr = ['T','S','rho','u','v']

for m in run_incr:
    for j in yr_incr:
        for k in mo_incr:
            for i in lat_incr:
                wed.transect([i,i],lonW,var_incr,j,k,varlim,run=m)
            for i in lon_incr:
                wed.transect(latS,[i,i],var_incr,j,k,varlim,run=m)
