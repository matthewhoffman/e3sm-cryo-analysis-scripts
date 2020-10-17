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
latS = [74,75]
lonW = [30,30]
startyr = 70
endyr = 101
yr = 93
m = 3

i=1
wed.profile(['T','S','rho'],'ISMF',yr,yr,latS[i],lonW[i],runcmp=True)#,mo=m)
wed.profile(['u','v'],'ISMF',yr,yr,latS[i],lonW[i],runcmp=True,mo=m)
