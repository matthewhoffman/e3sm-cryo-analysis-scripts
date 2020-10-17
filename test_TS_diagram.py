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
latS = 76 + 2.75/60
lonW = 30 + 59.65/60
startyr = 90#,60,70,80,90]
endyr = 80
dt = 9 
zeval = 20

#run_incr = 'ISMF'
#wed.TS_diagram(run_incr,latS,lonW,startyr,startyr+dt,z=zeval,zab=True,zall=True,seasonal=True)
#
#run_incr = 'ISMF-noEAIS'
#wed.TS_diagram(run_incr,latS,lonW,startyr,startyr+dt,z=zeval,zab=True,zall=True,seasonal=True)

run_incr = 'ISMF-3dGM'
wed.TS_diagram(run_incr,latS,lonW,startyr,startyr+dt,z=zeval,zab=True,zall=True,seasonal=True)

