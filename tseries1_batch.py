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

run_incr = ['ISMF']#,'ISMF-noEAIS']#,'ISMF-3dGM']
#var_incr = ['taux','tauy']
#var_incr = ['T','S','rho','u','v']
#var_incr = ['T','S']
#var_incr = ['unormal']
var_incr = ['rho']

#year_range = [70,101]
year_range = [90,91]

#print(len(var_incr))
#for i in run_incr:
wed.tseries1(run_incr,var_incr,year_range=year_range,
             #option = 'coord',placename = 'S4E',#'M31W',
             placename = 'gyre_interior',
             print_to_file=True,create_figure=True,
             ztop_pyc = [True], zbottom_pyc = [False],
             year_overlay=False)
#             zrange=[0,20],zab=True,
#        wed.tseries1(i,var_incr,latS,lonW,j,j,z=zeval,zab=True,runcmp=False,velocity_vector=True)
#        wed.tseries1(i,var_incr,latS,lonW,j,j+dt,z=zeval,zab=False,runcmp=True)

