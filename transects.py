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
lat_incr = np.arange(74,83,2)#1)
lon_incr = np.arange(30,40,2)#1)
varlim = True
yr_incr = np.arange(60,101,10)
#yr_incr = np.arange(64,91,2)
mo_incr = np.arange(1,13,3)
run_incr = ['ISMF','ISMF-noEAIS']#
var_incr = ['T','S','rho','u','v']

for m in run_incr:
    for j in yr_incr:
        #for k in mo_incr:
        wed.transect([76,76],[39,27],var_incr,j,11,varlim,run=m,new=True)
        #wed.transect([75,75],[36,25],var_incr,j,4,varlim,run=m)
        #wed.transect([78,78],[44,35],var_incr,j,1,varlim,run=m)
#for m in run_incr:
#    for j in yr_incr:
#        for k in mo_incr:
#            for i in lat_incr:
#                wed.transect([i,i],lonW,var_incr,j,k,varlim,run=m)
#            for i in lon_incr:
#                wed.transect(latS,[i,i],var_incr,j,k,varlim,run=m)
