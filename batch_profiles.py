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
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm
import cmocean
from math import pi
import weddell_mod as wed

run_incr = ['ISMF','ISMF-noEAIS']
latS = 76 + 2.75/60
lonW = 30 + 59.65/60
#latS = [74,75]
#lonW = [30,30]
startyr = 70
#endyr = 101
endyr = 71

for yr in np.arange(startyr,endyr+1,1):
    wed.profile(['T','S','rho'],'ISMF',yr,yr,latS,lonW,runcmp=True,runcmpname='ISMF-3dGM')
    wed.profile(['u','v'],'ISMF',yr,yr,latS,lonW,runcmp=True,runcmpname='ISMF-3dGM')
#    for i,l in enumerate(latS):
#        wed.profile(['T','S','rho'],'ISMF',yr,yr,latS[i],lonW[i],runcmp=True)
#        wed.profile(['u','v'],'ISMF',yr,yr,latS[i],lonW[i],runcmp=True)
#        for run in run_incr:
#            wed.profile(['T','S','rho'],run,yr,yr,latS[i],lonW[i])
#            wed.profile(['u','v'],run,yr,yr,latS[i],lonW[i])
