#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

from math import pi,sqrt,nan,atan,cos,floor,ceil
import matplotlib as pltlib
import numpy as np
import gsw
from data_config import *

pltlib.rc_file('rcparams.txt', use_default_template=True)

printformat = 'png'

m3ps_to_Sv = 1e-6 # m^3/sec flux to Sverdrups
kg_to_Gt = 1e-12
s_to_yr = 365*3600*24
bad_data = -1e33
deg2rad  = pi/180.
lat_N = -50 # northern limit of domain
rho_ice = 917
    
legloc = 'upper right'
#legloc = 'upper left'
bboxanchor = (0.0,-0.25)

#Ryan et al. 2017
season = ['summer','fall','winter','spring'] 
summer_color = '#DAF5D0'#(11,0,15,4)
fall_color   = '#F2D0DC'#(0,14,9,5)
winter_color = '#FFFFEB'#(0,0,8,0) 
spring_color = '#D1D5EB'#(11,9,0,8) 
season_color = [summer_color,fall_color,winter_color,spring_color]
t_season = [0,
            (2 + (15/31))/12,
            (6 + (5/30))/12,
            (9 + (5/31))/12]

set_dpi = 100

# This is what will be used for plotting
runtitle = runname

run_tipping_year = [71, nan, nan]
run_color = ['k',
             '#7570b3', # purple
             '#1b9e77'] # green


vartitle = ['z','z_pyc','T','S','rho','u','v','U','taux','tauy','tau','sea_ice_fw_flux','land_ice_fw_flux',
            'mean','ssh','ssh_cmp','curl','zice','u_barotropic','u_baroclinic','F_barotropic','F_baroclinic_pos']
varlabel = [str() for i in vartitle]
varlabel[vartitle.index('z')] = 'Depth (m)'
varlabel[vartitle.index('T')] = r'$T \: (^{\circ}C)$'
varlabel[vartitle.index('S')] = r'$S \: (PSU)$'
varlabel[vartitle.index('u')] = 'Velocity +East (m/s)'
varlabel[vartitle.index('v')] = 'Velocity +North (m/s)'
varlabel[vartitle.index('U')] = 'Velocity magnitude (m/s)'
varlabel[vartitle.index('rho')] = r'$\sigma_0 \: (kg \: m^{-3})$'
varlabel[vartitle.index('tau')] = 'Surface wind stress'
varlabel[vartitle.index('sea_ice_fw_flux')] = 'Sea ice freshwater flux ($kg s^{-1}$)'
varlabel[vartitle.index('land_ice_fw_flux')] = 'Land ice freshwater flux ($kg s^{-1}$)'
varlabel[vartitle.index('ssh')] = 'Sea Surface Height (m)'
varlabel[vartitle.index('curl')] = 'curl'
varlabel[vartitle.index('z_pyc')] = 'Pycnocline depth (m)'
varlabel[vartitle.index('mean')] = 'Pycnocline depth (m)'
varlabel[vartitle.index('F_barotropic')] = 'Barotropic Flux (Sv)'
varlabel[vartitle.index('F_baroclinic_pos')] = 'Baroclinic Flux (Sv)'
varlabel[vartitle.index('u_barotropic')] = 'u barotropic (m/s)'
varlabel[vartitle.index('u_baroclinic')] = 'u baroclinic (m/s)'
#varlabel[vartitle.index('u_barotropic')] = 'Velocity, barotropic (m/s)'
#varlabel[vartitle.index('u_baroclinic')] = 'Velocity, baroclinic (m/s)'
surfvar = ['taux','tauy','tau','ssh','sea_ice_fw_flux','land_ice_fw_flux']

# TODO fix assignment
varname  = ['depth','',
            'timeMonthly_avg_activeTracers_temperature',
            'timeMonthly_avg_activeTracers_salinity',
            'timeMonthly_avg_potentialDensity',
            'timeMonthly_avg_velocityZonal',
            'timeMonthly_avg_velocityMeridional',
            '','timeMonthly_avg_windStressZonal',
            'timeMonthly_avg_windStressMeridional',
            '','timeMonthly_avg_ssh','landIceDraft']
varname[vartitle.index('sea_ice_fw_flux')] = 'timeMonthly_avg_seaIceFreshWaterFlux'
varname[vartitle.index('land_ice_fw_flux')] = 'timeMonthly_avg_landIceFreshwaterFlux'
vartype = ['scalar' for i in vartitle]
vartype[vartitle.index('u')] = 'velocity'
vartype[vartitle.index('v')] = 'velocity'

# variable axis limits 
# choices for on the continental shelf
varmin = [str() for i in vartitle]
varmax = [str() for i in vartitle]
varmin[vartitle.index('z')] = -1800
varmax[vartitle.index('z')] = -100
varmin[vartitle.index('T')] = -2.1
varmax[vartitle.index('T')] = 1.7#-0.5
varmin[vartitle.index('S')] = 34.2#33.5
varmax[vartitle.index('S')] = 34.7#34.6
varmin[vartitle.index('rho')] = 27.40
varmax[vartitle.index('rho')] = 27.85
#varmin[vartitle.index('rho')] = 32.20
#varmax[vartitle.index('rho')] = 32.60
varmin[vartitle.index('u')] = -0.04
varmax[vartitle.index('u')] =  0.04
varmin[vartitle.index('v')] = -0.04
varmax[vartitle.index('v')] =  0.04
varmin[vartitle.index('U')] = -0.02
varmax[vartitle.index('U')] =  0.02
varmin[vartitle.index('tau')] = -0.003
varmax[vartitle.index('tau')] =  0.003
varmin[vartitle.index('taux')] = -0.003
varmax[vartitle.index('taux')] =  0.003
varmin[vartitle.index('tauy')] = -0.003
varmax[vartitle.index('tauy')] =  0.003
varmin[vartitle.index('ssh')] = -4
varmax[vartitle.index('ssh')] = -1
varmin[vartitle.index('ssh_cmp')] = -0.25
varmax[vartitle.index('ssh_cmp')] =  0.25
varmin[vartitle.index('u_barotropic')] = -0.0125
varmax[vartitle.index('u_barotropic')] = 0.0125
varmin[vartitle.index('u_baroclinic')] = -0.04
varmax[vartitle.index('u_baroclinic')] = 0.04
# choices for off the continental shelf
#varmin = [-4000,-2, 33.5, 1027.25, -0.03, -0.02, -0.03,-0.003,-0.003,-0.003,-5.]
#varmax = [-100, 2, 34.9, 1028   ,  0.03,  0.02, 0.03,0.003,0.003,0.003,0.]

# increment for variables
dvar = [str() for i in vartitle]
dvar[vartitle.index('T')] = 0.2
dvar[vartitle.index('S')] = 0.01
dvar[vartitle.index('rho')] = 0.01
dvar[vartitle.index('u')] = .001
dvar[vartitle.index('v')] = .001

# colormap corresponding to var
varcmap = ['cmo.speed' for i in vartitle]
varcmap[vartitle.index('z')] = 'cmo.deep'
varcmap[vartitle.index('z_pyc')] = 'cmo.deep'
varcmap[vartitle.index('zice')] = 'cmo.deep'
varcmap[vartitle.index('T')] = 'cmo.thermal'
varcmap[vartitle.index('S')] = 'cmo.haline'
varcmap[vartitle.index('rho')] = 'cmo.dense'
varcmap[vartitle.index('u')] = 'cmo.balance'
varcmap[vartitle.index('v')] = 'cmo.balance'
varcmap[vartitle.index('u_barotropic')] = 'cmo.balance'
varcmap[vartitle.index('u_baroclinic')] = 'cmo.balance'
varcmap[vartitle.index('tau')] = 'cmo.curl'
varcmap[vartitle.index('taux')] = 'cmo.balance'
varcmap[vartitle.index('tauy')] = 'cmo.balance'
varcmap[vartitle.index('U')] = 'cmo.speed'
varcmap[vartitle.index('ssh')] = 'cmo.speed'

# fontsize
fs=14
# linewidth
lw1 = 1.0

TSpolygon_Hattermann2018 = [[34.45,-0.5],[34.5,0],[34.6,-1],[34.55,-1.5]]
#TSpolygon_Hattermann2018_edit = [[34.20,-0.75],[34.25,0.25],[34.35,-0.75],[34.30,-1.75]]
TSpolygon_Hattermann2018_edit = [[34.10,0.25],[34.15,1.25],[34.35,-0.75],[34.30,-1.75]]
TSpolygon_ISMF = [[34.45,0.5],[34.5,0],[34.5,-1],[34.0,-1.5]]

#TODO
TSpolygon_CDW  = [[34.0, 0.0], [40.0, 0.0], [40.0, 5.0], [34.0, 5.0]] 
TSpolygon_mCDW = [[34.0, -1.5], [40.0, -1.5], [40.0, 0.0], [34.0, 0.0]] 

m_Tfreezing = -0.0575
b_Tfreezing = 0.0901
# Tfpreezing = m_Tfreezing * S + b_Tfreezing
TSpolygon_AASW = [[20.0, 5.0], [20.0, m_Tfreezing * 20.0 + b_Tfreezing],
                  [34.0, m_Tfreezing * 34.0 + b_Tfreezing], [34.0, 5.0]]
TSpolygon_LSSW = [[34.0, -1.5], [34.0, m_Tfreezing * 34.0 + b_Tfreezing],
                  [34.5, m_Tfreezing * 34.5 + b_Tfreezing], [34.5, -1.5]]
TSpolygon_HSSW = [[34.5, -1.5], [34.5, m_Tfreezing * 34.5 + b_Tfreezing],
                  [40.0, m_Tfreezing * 40.0 + b_Tfreezing], [40.0, -1.5]]
TSpolygon_ISW  = [[30.0, -3.0], [30.0, m_Tfreezing * 30.0 + b_Tfreezing],
                  [40.0, m_Tfreezing * 40.0 + b_Tfreezing], [40.0, -3.0]]
