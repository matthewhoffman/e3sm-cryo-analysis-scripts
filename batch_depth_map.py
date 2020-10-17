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
var = ['T']
c = 0.8
#var = ['rho']
#c = 1027.6
var = ['z_pyc']
run = ['ISMF']#,'ISMF-noEAIS']
yr = [70]
mo = [1]
#yr = np.arange(90,95,1)
#mo = np.arange(1,13,1)
#yr = np.arange(70,101,1)
#mo = np.arange(1,13,1)
#loc = 'fris'
loc = 'wed_pyc_Ryan'
#loc = 'frisEAcoast'

#wed.plot_mesh_var('zice')
for k in yr:
    for m in mo:
        for j in var:
            #wed.plot_surf_var(j,k,m,z=zi,run=run[0],zab=False,level=c,runcmp=True,locname=loc,varlim=True)
            for i in run:
               wed.plot_surf_var(j,k,m,z=zi,run=[i],zab=False,
                                 locname=loc,varlim=False,overwrite=True)
               #wed.plot_surf_var(j,k,m,z=0,run=i,zab=False,runcmp=False,locname=loc,varlim=True)
               #wed.plot_surf_var(j,k,m,z=0,run=i,zab=True,runcmp=False,locname=loc,varlim=True)

