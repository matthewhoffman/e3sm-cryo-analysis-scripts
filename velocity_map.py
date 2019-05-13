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
from math import pi
from matplotlib.colors import LogNorm

#def plot_surf_var(var_incr,z,yr,mo,run='ISMF-noEAIS'):# user-defined parameters
var = 'ssh'
run = 'ISMF'
vartitle = ['T','S','rho','u','v','U','tau','ssh']
varname  = ['timeMonthly_avg_activeTracers_temperature',
                'timeMonthly_avg_activeTracers_salinity',
                'timeMonthly_avg_potentialDensity',
                'timeMonthly_avg_velocityZonal',
                'timeMonthly_avg_velocityMeridional',
                '','',
                'timeMonthly_avg_ssh']
varmin = [-2.5,33.5,1027,-0.01,-0.01]
varmax = [1,34.4,1027.75,0.01,0.01]
varcmap = ["viridis","viridis","viridis","coolwarm","coolwarm"]

yr = 70
mo = 1

lat_N = -50 # northern limit of domain

run = 'ISMF'
if run == 'ISMF':
    path = '/usr/projects/climate/mhoffman/e3sm/cryo-campaign-data'
#    path = '/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/run'
elif run == 'ISMF-noEAIS':
    path = '/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run'
else:
    print('run does not exist')

#fmesh = netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
fmesh = netCDF4.Dataset('../oEC60to30v3wLI60lev.171031.nc')
filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.{1:04d}-{2:02d}-01.nc'.format(path, yr, mo)
f = netCDF4.Dataset(filename, 'r')
#if cmp:

# constants
deg2rad  = pi/180.
bad_data = -1e33
bad_data2 = 0.

# import variables from file
latCell  = fmesh.variables['latCell'][:]
lonCell  = fmesh.variables['lonCell'][:]
idxCell  = fmesh.variables['indexToCellID'][:]
xCell    = fmesh.variables['xCell'][:]
yCell    = fmesh.variables['yCell'][:]
depths   = fmesh.variables['refBottomDepth'][:]
kmax     = fmesh.variables['maxLevelCell'][:]
zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
zh       = fmesh.variables['layerThickness'][0,:]
icemask  = fmesh.variables['landIceMask'][:]
zice     = fmesh.variables['landIceDraft'][0,:]
ssh      = f.variables['timeMonthly_avg_ssh'][0,:]

# northern limit for subplots
logical_N = (xCell<2.5e6) & (xCell > 0) & (yCell < 0) & (yCell > -1.5e6) & (latCell < lat_N*deg2rad) 

zbottom  = np.zeros(zh.shape)
zbottom[:,-1] = zmax
for i in range(len(depths)-2,-1,-1):
    zbottom[:,i] = zbottom[:,i+1] + zh[:,i+1]
zmid = zbottom + np.multiply(0.5,zh)
zssh = zbottom[:,0] + zh[:,0]

if var == 'U':
    u = f.variables['timeMonthly_avg_velocityZonal'][0,logical_N,:]
    v = f.variables['timeMonthly_avg_velocityMeridional'][0,logical_N,:]
elif var == 'tau':
    u = f.variables['timeMonthly_avg_windStressZonal'][0,logical_N]
    v = f.variables['timeMonthly_avg_windStressMeridional'][0,logical_N]
elif var == 'ssh':
    data = ssh
else:
    data = f.variables[varname[vartitle.index(var)]][0,logical_N,:]

# define plot location
if (var == 'U'):
    for i in range(0,len(u),1):
        uz[i] = data[i,np.argmin(np.abs(np.subtract(zmid[i,:],z)))]
        vz[i] = data[i,np.argmin(np.abs(np.subtract(zmid[i,:],z)))]
        dataz = np.sqrt(np.add(np.square(uz),np.square(vz)))
        heading = np.divide(np.arctan2(uz,vz),deg2rad)
elif (var == 'tau'):
        dataz = np.sqrt(np.add(np.square(u),np.square(v)))
        heading = np.divide(np.arctan2(u,v),deg2rad)
else:
    dataz = np.zeros(len(data),)
    for i in range(0,len(data),1):
        dataz[i] = data[i,np.argmin(np.abs(np.subtract(zmid[i,:],z)))]


# plots
size = 30#size = 60
fig = plt.figure()
#cntr1 = ax1.tricontourf(yCell[logical_N], xCell[logical_N], 
#                        U, 20, cmap="viridis")
cntr1 = plt.scatter(yCell[logical_N], xCell[logical_N], 
                    s=size, c=dataz, marker = '.',cmap="viridis",
                    vmin = 0.01, vmax = 0.05, norm = LogNorm())
ax = fig.gca()
ax.set_aspect('equal')
gl1 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                       zice[logical_N].flatten(), [-1e-10], colors = 'k', linewidths = 2)
cbar1 = fig.colorbar(cntr1)
cbar1.set_label(var)
fig.tight_layout()

if (var == 'U') | (var == 'tau'):
    fig = plt.figure()
    #cntr2 = ax2.tricontourf(yCell[logical_N], xCell[logical_N], 
    #                        heading, 20, cmap="twilight", vmin = -180, vmax = 180)
    cntr2 = plt.scatter(yCell[logical_N], xCell[logical_N], 
                        s=size, c=heading,  marker = '.',cmap="twilight", vmin = -180, vmax = 180)
    gl2 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                           zice[logical_N].flatten(), [-1e-10], colors = 'k', linewidths = 2)
    ax = fig.gca()
    ax.set_aspect('equal')
    cbar2 = fig.colorbar(cntr2)
    cbar2.set_label('Degrees from N')    
    fig.tight_layout()
    plt.savefig('velocity_heading' + str(yr) + '-' + str(mo) + '.png')
#plt.close()
