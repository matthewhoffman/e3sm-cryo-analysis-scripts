#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.signal
import cmocean
from matplotlib import cm
from math import pi
from matplotlib.colors import LogNorm
import weddell_mod as wed
#cmaps = cmocean.cm.cmap_d

# PARAMETERS

zi = 0.;
#zi = 100.;
var = ['U']
#var = ['tau']
#var = ['ssh']
#var = ['rho']
#run = ['ISMF','ISMF-noEAIS']
run = ['ISMF','ISMF-3dGM']
yr = np.arange(90,101,1)
mo = np.arange(1,13,1)
#yr = np.arange(90,91,1)
#mo = np.arange(1,2,1)
plottype='quiver'
#loc = 'fris'
loc = 'frisEAcoast'

for k in yr:
    for m in mo:
        for j in var:
            for i in run:
                wed.plot_surf_var(j,k,m,z=zi,run=[i],zab=False,plottype=plottype,
                              locname=loc,varlim=True,overwrite=True)
            wed.plot_surf_var(j,k,m,z=zi,run=run,zab=False,plottype=plottype,
                              locname=loc,varlim=True,overwrite=True)

