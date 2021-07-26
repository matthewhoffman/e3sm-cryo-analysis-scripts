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

run_incr = ['ISMF','ISMF-noEAIS','ISMF-3dGM','ISMF-noDIB']
#lat = -70
#lon = 360-20
year_range = [70,70]
varlist = ['T','S','rho']

for i in range(20,180,10):
    wed.profile(run_incr,varlist,[i,i],mo=8,placename='N31W',maxDepth=-800)
    wed.profile(run_incr,varlist,[i,i],mo=2,placename='N31W',maxDepth=-800)
#    wed.profile(['u','v'],'ISMF',yr,yr,latS,lonW,runcmp=True,runcmpname='ISMF-3dGM')
#    for i,l in enumerate(latS):
#        wed.profile(['T','S','rho'],'ISMF',yr,yr,latS[i],lonW[i],runcmp=True)
#        wed.profile(['u','v'],'ISMF',yr,yr,latS[i],lonW[i],runcmp=True)
#        for run in run_incr:
#            wed.profile(['T','S','rho'],run,yr,yr,latS[i],lonW[i])
#            wed.profile(['u','v'],run,yr,yr,latS[i],lonW[i])
