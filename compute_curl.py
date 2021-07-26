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

#run_incr = ['ISMF','ISMF-noEAIS']
#run_incr = ['ISMF','B-ISMF']
run_incr = ['ISMF-noDIB','ISMF-3dGM']
startyr = 10
endyr = 180
#endyr = 21
wed.calc_stresscurl_t(run_list=run_incr,year_range=[startyr,endyr],region='wedwang')
#wed.plot_stresscurl_t('windstresscurl_wedwang_50-101',run_incr,year_range=[70,81])
#wed.plot_stresscurl_t('ISMF_ISMF-noEAIS_windstresscurl_wed_pyc_Ryan_10-180',run_incr,year_range=[20,170],region='wed_pyc_Ryan')
#wed.plot_stresscurl_t('ISMF_B-ISMF_windstresscurl_wedwang_20-100',run_incr,year_range=[20,100],region='wedwang',plot_difference = True)
#wed.fluxgate('trough_ice',['T'],50,1,overwrite=True)
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [70.,80.])
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [90.,100.])
