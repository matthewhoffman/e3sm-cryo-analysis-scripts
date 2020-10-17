#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

from math import pi,sqrt,nan,atan,cos,floor,ceil
import matplotlib as pltlib
import numpy as np

pltlib.rc_file('rcparams.txt', use_default_template=True)

printformat = 'png'

m3ps_to_Sv = 1e-6 # m^3/sec flux to Sverdrups
bad_data = -1e33
deg2rad  = pi/180.
lat_N = -50 # northern limit of domain
    
cells_trough_ice_lat = [33000,77033,11788,129098,49348] # index is value - 1, located just under ice
cells_trough_shelf_lat = [165785,166563,130260,191987,85569] # index is value - 1, located just under ice
#cells_trough_shelf_lat = [224105,165785,166563,130260,191987] # index is value - 1, located just under ice
cells_trough_shelf2_lat = [152383,144943,67352] # index is value - 1, located just under ice

latmin = -80
latmax = -80
lonmin = 360-48
lonmax = 360

region_name = ['frisEAcoast','fris','gyre_interior','trough_shelf','trough_ice','M31W',
               'S4E','wedwang','wed_pyc_Ryan'] # formerly loc
# M31W and S4E from Ryan et al. 2017 JGR Oceans
region_title = ['' for i in region_name] # loc_name
#TODO fix allocation
region_coordbounds = np.zeros((len(region_name),2,2)) 
region_xybounds    = np.zeros((len(region_name),2,2)) 
region_zbounds     = np.zeros((len(region_name),2)) 
for i in region_name:
    region_zbounds[region_name.index(i)] = [-9999.,20.]

region_title[region_name.index('gyre_interior')] = 'off-shelf'
region_title[region_name.index('trough_shelf')] = 'trough, shelf break'
region_title[region_name.index('trough_shelf')] = 'trough, shelf break'
region_title[region_name.index('trough_ice')] = 'trough, under ice shelf front'
region_title[region_name.index('M31W')] = 'M31W'
region_title[region_name.index('S4E')] = 'S4E'

region_coordbounds[region_name.index('M31W')] = ([360-30.994,360-30.994],
                                                 [-76.046,-76.046])
region_coordbounds[region_name.index('S4E')] = ([360-30.58,360-30.58],
                                               [-74.62,-74.62])
region_xybounds   [region_name.index('frisEAcoast')]  = [0e6,2e6],[-2e6,0] 
region_coordbounds[region_name.index('frisEAcoast')]  = [300,360],[-180,-60]
region_xybounds   [region_name.index('fris')]         = [0.2e6,1.5e6],[-1.5e6,-0.4e6] 
region_coordbounds[region_name.index('fris')]         = [-1,361],[-181,181]
region_xybounds   [region_name.index('wedwang')]      = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wedwang')]      = [305,315],[-74,-65] 
region_xybounds   [region_name.index('wed_pyc_Ryan')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wed_pyc_Ryan')] = [330,360],[-75,-69] 
region_xybounds   [region_name.index('trough_shelf')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('trough_shelf')] = [325,330],[-74.75,-74.75]
region_xybounds   [region_name.index('trough_ice')]   = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('trough_ice')]   = [314,325],[-78.1,-78.1]
region_xybounds   [region_name.index('gyre_interior')]= region_xybounds[region_name.index('wed_pyc_Ryan')]
region_coordbounds[region_name.index('gyre_interior')]= region_coordbounds[region_name.index('wed_pyc_Ryan')]
region_zbounds    [region_name.index('gyre_interior')]= [-4500,-2500]

loc_ptsize = [40 for i in region_name]

legloc = 'upper left'
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

runtitle = ['ISMF-control',
           'ISMF-noEAIS',
           'ISMF-3dGM']
runname = ['ISMF',
           'ISMF-noEAIS',
           'ISMF-3dGM']
run_color = ['cornflowerblue',
             'k',
             'darkorange']
savepath_nersc = '/global/homes/c/cbegeman/weddell_output/'
runpath = [
#          '/global/cfs/cdirs/m3412/simulations/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
           '/global/cfs/cdirs/m3412/simulations/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
#           '/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
           
#            '/global/cfs/cdirs/m3412/simulations/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/archive/ocn/hist',
             '/global/cfs/cdirs/m3412/simulations/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist/',
           #'/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run',
#           '/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist',

           '/global/cscratch1/sd/sprice/acme_scratch/cori-knl/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/archive/ocn/hist']
meshpath = ['/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc',
            '/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc',
            '/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc']

vartitle = ['z','z_pyc','T','S','rho','u','v','U','taux','tauy','tau',
            'ssh','ssh_cmp','curl','zice','u_barotropic','u_baroclinic']
varlabel = [str() for i in vartitle]
varlabel[0:7] = ['Depth (m)','T (degC)','S (PSU)','Pot. Density (kg/m3)',
                 'Velocity +East (m/s)','Velocity +North (m/s)',
                 'Velocity magnitude (m/s)']
varlabel[vartitle.index('T')] = 'T'
varlabel[vartitle.index('S')] = 'S'
varlabel[vartitle.index('rho')] = 'Pot density'
varlabel[vartitle.index('tau')] = 'Surface stress'
varlabel[vartitle.index('ssh')] = 'Sea Surface Height (m)'
varlabel[vartitle.index('curl')] = 'curl'
varlabel[vartitle.index('z_pyc')] = 'Pycnocline depth (m)'
varlabel[vartitle.index('u_barotropic')] = 'u barotropic (m/s)'
varlabel[vartitle.index('u_baroclinic')] = 'u baroclinic (m/s)'
#varlabel[vartitle.index('u_barotropic')] = 'Velocity, barotropic (m/s)'
#varlabel[vartitle.index('u_baroclinic')] = 'Velocity, baroclinic (m/s)'
surfvar = ['taux','tauy','tau','ssh']

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

vartype = ['scalar' for i in vartitle]
vartype[vartitle.index('u')] = 'velocity'
vartype[vartitle.index('v')] = 'velocity'

# variable axis limits 
# choices for on the continental shelf
varmin = [str() for i in vartitle]
varmax = [str() for i in vartitle]
varmin[vartitle.index('z')] = -1800
varmax[vartitle.index('z')] = -100
varmin[vartitle.index('T')] = -2.0
varmax[vartitle.index('T')] = 1.0#-0.5
varmin[vartitle.index('S')] = 33.5#34.2
varmax[vartitle.index('S')] = 34.6
varmin[vartitle.index('rho')] = 1027.20#.40
varmax[vartitle.index('rho')] = 1027.75
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
dvar[vartitle.index('T')] = 0.1
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
fs=16
# linewidth
lw1 = 1

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
