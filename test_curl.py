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
startyr = 50
endyr = 101
yr = 93
m = 3
i=1
teval = [80.,90.]
locname = 'trough_shelf' 
#wed.pick_transect(transect_name=locname,overwrite=True,option='by_index',scope_name='fris')
wed.plot_stresscurl_t('windstresscurl_wedwang_50-101',run_incr,year_range=[70,81])
#wed.fluxgate('trough_ice',['T'],50,1,overwrite=True)
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [70.,80.])
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [90.,100.])
#wed.fluxgate(locname,yrrange=[startyr,endyr],overwrite=False,runcmp=True)#,plotting=True)
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr))
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = teval)

#locname = 'trough_ice' 
#wed.pick_transect(transect_name=locname,overwrite=True,option='by_index',scope_name='fris')
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = teval)
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = [90.,100.])
