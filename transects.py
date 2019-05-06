#!/usr/bin/env python

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm
from math import pi

# user-defined parameters
var = 'T'
dir_transect = 'y' # defines whether profile should be sorted by x or y
#lat_transect = [-75, -75]
#lon_transect = [360-70, 360]

lat_transect = [-78, -78]
lon_transect = [360-70, 360]

#lat_transect = [-85, -75]
#lon_transect = [360-60, 360-60]

#lon_transect = [360-40, 360-40]
lat_N = -75 # northern limit of domain

fmesh = netCDF4.Dataset('../oEC60to30v3wLI60lev.171031.nc')
#fmesh=netCDF4.Dataset('/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc')
f = netCDF4.Dataset('/usr/projects/climate/mhoffman/e3sm/cryo-campaign-data/mpaso.hist.am.timeSeriesStatsMonthly.0070-01-01.nc')

# constants
deg2rad  = pi/180.
vartitle = ['T','S','rho']
varname  = ['timeMonthly_avg_activeTracers_temperature',
            'timeMonthly_avg_activeTracers_salinity',
            'timeMonthly_avg_potentialDensity']
dlat = 0.15 # at 30km resolution, distance between cells in latitude space
dlon = .75
bad_data = -1e33

# inport variables from file
latCell  = fmesh.variables['latCell'][:]
lonCell  = fmesh.variables['lonCell'][:]
idxCell  = fmesh.variables['indexToCellID'][:]
xCell    = fmesh.variables['xCell'][:]
yCell    = fmesh.variables['yCell'][:]
depths   = fmesh.variables['refBottomDepth'][:]
zmax     = fmesh.variables['bottomDepth'][:]
icemask  = fmesh.variables['landIceMask'][:]
zice     = fmesh.variables['landIceDraft'][0,:]
data     = f.variables[varname[vartitle == var]][:]
ssh      = f.variables['timeMonthly_avg_ssh'][0,:]

# calculate z from depths
z        = np.zeros(depths.shape)
z[0]     = -0.5 * depths[0]
z[1:]    = -0.5 * (depths[0:-1] + depths[1:])

# define line of constant latitude
if lat_transect[0] == lat_transect[1]:
    lat_transect[0] = lat_transect[0] - dlat
    lat_transect[1] = lat_transect[1] + dlat
if lon_transect[0] == lon_transect[1]:
    lon_transect[0] = lon_transect[0] - dlon
    lon_transect[1] = lon_transect[1] + dlon

# define plot location

# northern limit for subplots
logical_N = (latCell < lat_N*deg2rad) & (xCell > 0)

# indices of transect
logical_trans = ( (latCell > lat_transect[0]*deg2rad) & 
                  (latCell < lat_transect[1]*deg2rad) &
                  (lonCell > lon_transect[0]*deg2rad) & 
                  (lonCell < lon_transect[1]*deg2rad)   )
idx_trans = np.where(logical_trans)[0]

y1 = np.min(yCell[logical_trans])
x1 = np.max(xCell[logical_trans])
temp = np.sqrt( np.square(yCell[logical_trans] - y1) + 
                np.square(xCell[logical_trans] - x1)   )
idxsort_trans = idx_trans[temp.argsort()]
ysort_trans = yCell[idxsort_trans]
xsort_trans = xCell[idxsort_trans]
dist_trans = temp[temp.argsort()]

# distance along transect

# create mesh variables for plotting
ymesh,zmesh = np.meshgrid(dist_trans, z);
sshmesh,temp = np.meshgrid(ssh[idxsort_trans],z)
data_trans = np.transpose(data[0,idxsort_trans,:])

# Ideally, mask both bad data and ice shelf
#data_trans_masked = np.ma.masked_where( (data_trans < bad_data) |
#                                        (zmesh > sshmesh), data_trans)
data_trans_masked = np.ma.masked_where( (data_trans < bad_data), data_trans)

zadj = np.add(zmesh,zice[logical_trans])
# plots

# show profile line across cells
fig = plt.figure()
plt.plot(yCell[logical_N],     xCell[logical_N], 'k.')
plt.plot(yCell[logical_trans], xCell[logical_trans], 'r.')
plt.axis('equal')
plt.savefig('test_domain.png')
#plt.close()

fig = plt.figure()
plt.pcolormesh(np.divide(ymesh,1e3), zmesh, data_trans_masked)#, vmin=-2.1, vmax=1.6, cmap = 'viridis')
plt.plot(np.divide(dist_trans,1e3), zice[logical_trans], color = 'red', marker = '.', linestyle = '-')
plt.plot(np.divide(dist_trans,1e3), -1.*zmax[logical_trans], color = 'black', marker = '.', linestyle = '-')
plt.ylim((min(-1.*zmax[logical_trans]),0))
plt.colorbar()
plt.xlabel('Distance (km)')
plt.ylabel('Depth (m)')
plt.title(var)
plt.savefig('test_transect.png')

fig = plt.figure()
plt.scatter(yCell[logical_N], xCell[logical_N], s=5, c = zice[logical_N], cmap = 'viridis') # c=lonCell[idx1]/deg2rad
plt.colorbar()
plt.axis('equal')
plt.title('Land Ice Draft (m)')

fig = plt.figure()
plt.scatter(yCell[logical_N], xCell[logical_N], s=5, c = ssh[logical_N], cmap = 'viridis') # c=lonCell[idx1]/deg2rad
plt.colorbar()
plt.axis('equal')
plt.title('Average SSH (m)')

fig = plt.figure()
plt.scatter(yCell[logical_N], xCell[logical_N], s=5, c = zmax[logical_N], cmap = 'viridis') # c=lonCell[idx1]/deg2rad
plt.colorbar()
plt.axis('equal')
plt.title('Bottom depth (m)')

fig = plt.figure()
plt.scatter(yCell[logical_N], xCell[logical_N], s=5, c = ssh[logical_N] - zmax[logical_N], cmap = 'viridis') # c=lonCell[idx1]/deg2rad
plt.colorbar()
plt.axis('equal')
plt.title('Water column thickness (m)')
