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

run_incr = ['ISMF','ISMF-noEAIS','ISMF-3dGM']
latS = [74,75]
lonW = [30,30]
startyr = 70
endyr = 100 
i=1

#wed.pick_transect(transect_name=locname,overwrite=True,option='by_index',scope_name='fris')

#wed.plot_stresscurl_t('ISMF',yrrange=[startyr,endyr],overwrite=True,runcmp=True)
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [70.,80.])
#wed.plot_stresscurl_t_diff('windstresscurl_wedwang'+'_'+str(startyr)+'-'+str(endyr),tlim = [90.,100.])

#for r in run_incr:
#    wed.fluxgate(locname,yrrange=[startyr,endyr],morange = [0,12+1],overwrite=False,run=r,
#                 savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
#             runcmp=True,runcmpname='ISMF-3dGM')#,plotting=True)

#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr))

locname = 'trough_shelf' 
wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [70,80],
                    savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [90,100],
                    savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
locname = 'trough_ice' 
wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [70,80],
                    savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
wed.plot_fluxgate_t('ISMF_transect_u_'+locname+'_'+str(startyr)+'-'+str(endyr),
                    run_incr = run_incr, mode = 'barotropic-baroclinic',tlim = [90,100],
                    savepath = '/global/homes/c/cbegeman/weddell_output/fluxgate/')
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = teval)

#wed.pick_transect(transect_name=locname,overwrite=True,option='by_index',scope_name='fris')
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = teval)
#wed.plot_fluxgate_t('transect_flux_'+locname+'_'+str(startyr)+'-'+str(endyr),tlim = [90.,100.],runcmpname='ISMF-3dGM')
