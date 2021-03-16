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

#var_incr = ['taux','tauy']
#var_incr = ['T','S','rho','u','v']
var_incr = ['T','S','u','v']
#lat = -1*(76 + 2.75/60)
#lon = 360 - (30 + 59.65/60)
# for abysal profiles
lat = -73
lon = 360-30

year_range = [80,80]
dt = 1 
zeval = 20

#run_incr = 'ISMF'
run_incr = ['ISMF','ISMF-noEAIS']
#wed.TS_diagram(run_incr,lat,lon,startyr,startyr+dt,z=zeval,zab=False,zall=True,seasonal=True)
wed.TS_diagram(run_incr,year_range,
               zall=True,
               seasonal=False,plot_lines=False)
#wed.TS_diagram(lat,lon,startyr,startyr+dt,
#               run_list=run_incr,zall=True,
#               seasonal=False,plot_lines=False)
#
#run_incr = 'ISMF-noEAIS'
#wed.TS_diagram(run_incr,latS,lonW,startyr,startyr+dt,z=zeval,zab=True,zall=True,seasonal=True)

#run_incr = 'ISMF-3dGM'
#wed.TS_diagram(run_incr,latS,lonW,startyr,startyr+dt,z=zeval,zab=True,zall=True,seasonal=True)

