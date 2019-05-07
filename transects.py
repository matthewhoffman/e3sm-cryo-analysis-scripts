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

var = 'T'
lat_incr = np.arange(74,82.5,0.5)
lon_incr = np.arange(30,2,70)
lonW = [70, 0]
latS = [82, 74]
yr = 70
mo = 1
varlim = True
for i in lat_incr:
    wed.transect([i,i],lonW,'T',yr,mo,varlim)
    wed.transect([i,i],lonW,'S',yr,mo,varlim)
#    wed.transect([i,i],lonW,'rho',yr,mo,varlim)
#    wed.transect([i,i],lonW,'u',yr,mo,varlim)
#    wed.transect([i,i],lonW,'v',yr,mo,varlim)
for i in lon_incr:
    wed.transect(latS,[i,i],'T',yr,mo,varlim)
    wed.transect(latS,[i,i],'S',yr,mo,varlim)
#    wed.transect(latS,[i,i],'rho',yr,mo,varlim)
#    wed.transect(latS,[i,i],'u',yr,mo,varlim)
#    wed.transect(latS,[i,i],'v',yr,mo,varlim)
