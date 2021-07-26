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
from pick_from_mesh import pick_transect 

run_incr = ['ISMF','ISMF-noDIB','ISMF-noEAIS','ISMF-3dGM']
#latS = [74,75]
#lonW = [30,30]
startyr = 10
endyr = 180 
#i=1
savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/'
locname = 'wedwang_c'
#filename = 'ISMF_transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr)
#pick_transect(transect_name = locname,overwrite=True)
wed.fluxgate(locname,yrrange=[startyr,endyr],morange = [1,12+1],run_incr=run_incr,
             mode = 'barotropic-baroclinic', #mode = 'pos-neg',
             plot_map = False, plot_transect = False,
             overwrite=True,savepath = savepath)
#wed.plot_fluxgate_t(filename,tlim=[70,100],
#                    run_incr = run_incr, mode = 'pos-neg',
#                    var_incr = ['total'],savepath = savepath)

#wed.plot_stresscurl_t('ISMF',yrrange=[startyr,endyr],overwrite=True,runcmp=True)
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [70.,80.])
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [90.,100.])


#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr))

#locname = 'trough_shelf' 
#wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
#                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [70,80],
#                    savepath = savepath)
#wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
#                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [90,100],
#                    savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
#locname = 'trough_ice' 
#wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
#                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [70,80],
#                    savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
#wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
#                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [90,100],
#                    savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = teval)

#wed.pick_transect(transect_name=locname,overwrite=True,option='by_index',scope_name='fris')
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = teval)
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = [90.,100.],runcmpname='ISMF-3dGM')
