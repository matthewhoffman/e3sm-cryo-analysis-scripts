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

latS = [76,60]
lonW = [70,0]
#lat_incr = np.arange(74,83,2)#1)
#lon_incr = np.arange(30,40,2)#1)
varlim = True
yr_incr = np.arange(80,101,1)
mo_incr = np.arange(1,13,1)
run_incr = ['ISMF','ISMF-noEAIS']
var_incr = ['T','S','rho','u','v']

yr = 70
mo = 1
#wed.transect(latS,[30,30],['T'],yr,mo,varlim,run='ISMF',zlim='log')
#wed.transect(latS,[30,30],['T'],yr,mo,varlim,run='ISMF-noEAIS',zlim='log')
wed.transect(latS,[30,30],['T'],yr,mo,varlim,run='ISMF',zlim='log',runcmp=True)
