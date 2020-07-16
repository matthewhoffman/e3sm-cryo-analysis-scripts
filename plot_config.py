#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

from math import pi,sqrt,nan,atan,cos,floor,ceil
import matplotlib as pltlib

pltlib.rc_file('rcparams.txt', use_default_template=True)

m3ps_to_Sv = 1e-6 # m^3/sec flux to Sverdrups
bad_data = -1e33
#bad_data2 = 0.
deg2rad  = pi/180.
lat_N = -50 # northern limit of domain
    
cells_trough_ice_lat = [33000,77033,11788,129098,49348] # index is value - 1, located just under ice
cells_trough_shelf_lat = [165785,166563,130260,191987,85569] # index is value - 1, located just under ice
#cells_trough_shelf_lat = [224105,165785,166563,130260,191987] # index is value - 1, located just under ice
cells_trough_shelf2_lat = [152383,144943,67352] # index is value - 1, located just under ice

# location M31W from Ryan et al. 2017 JGR Oceans
latS = 76 + 2.75/60
lonW = 30 + 59.65/60

latmin = -80
latmax = -80
lonmin = 360-48
lonmax = 360

loc = ['frisEAcoast','fris','wedwang','trough_shelf','trough_ice','76S31W','75S31W']
loctitle = ['' for i in loc]
loctitle[loc.index('trough_shelf')] = 'trough, shelf break'
loctitle[loc.index('trough_ice')] = 'trough, under ice shelf front'
loctitle[loc.index('76S31W')] = 'M31W'
loctitle[loc.index('75S31W')] = 'S4E'

loc_xmin = [0e6,0.2e6,360-55,360-35,360-46]
loc_xmax = [2e6,1.5e6,360-45,360-30,360-35]
loc_ymin = [-2e6,-1.5e6,-74,-74.75,-78.1]
loc_ymax = [0,-0.4e6,-65,-74.75,-78.1]
loc_ptsize = [40,80]

legloc = 'upper left'
bboxanchor = (0.0,-0.25)

season = ['summer','fall','winter','spring'] 
summer_color = '#DAF5D0'#(11,0,15,4)
fall_color = '#F2D0DC'#(0,14,9,5)
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
runpath = ['/global/cfs/cdirs/m3412/simulations/20190225.GMPAS-DIB-IAF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
#/global/cfs/cdirs/m3412/simulations/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
#           '/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
           '/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist',
           '/global/cscratch1/sd/sprice/acme_scratch/cori-knl/20190819.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl.testNewGM/archive/ocn/hist']
           #'/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run']
meshpath = ['/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc',
            '/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc',
            '/project/projectdirs/e3sm/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc']

vartitle = ['z','T','S','rho','u','v','U','taux','tauy','tau',
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
varlabel[vartitle.index('u_barotropic')] = 'u barotropic (m/s)'
varlabel[vartitle.index('u_baroclinic')] = 'u baroclinic (m/s)'
#varlabel[vartitle.index('u_barotropic')] = 'Velocity, barotropic (m/s)'
#varlabel[vartitle.index('u_baroclinic')] = 'Velocity, baroclinic (m/s)'
surfvar = ['taux','tauy','tau','ssh']

# TODO fix assignment
varname  = ['depth','timeMonthly_avg_activeTracers_temperature',
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

TSpolygon_Hattermann2018 = [[34.45,0.5],[34.5,0],[34.6,-1],[34.55,-1.5]]
TSpolygon_ISMF = [[34.45,0.5],[34.5,0],[34.5,-1],[34.0,-1.5]]
