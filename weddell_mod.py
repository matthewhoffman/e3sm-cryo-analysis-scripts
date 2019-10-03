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
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from matplotlib import ticker

global bad_data, bad_data2, deg2rad, lat_N, runname, runpath, meshpath, vartitle, varname, varmin, varmax, varcmap, surfvar, dvar

bad_data = -1e33
bad_data2 = 0.
deg2rad  = pi/180.
lat_N = -50 # northern limit of domain
    
loc = ['frisEAcoast','fris']
loc_xmax = [2.5e6,1.5e6]
loc_ymin = [-1.5e6,-1.5e6]
loc_ptsize = [40,80]

set_dpi = 300.
    
runname = ['ISMF',
           'ISMF-noEAIS']
runpath = ['/global/cscratch1/sd/dcomeau/acme_scratch/cori-knl/20190225.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
           '/global/cscratch1/sd/hoffman2/acme_scratch/edison/archive/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/ocn/hist']
           #'/global/cscratch1/sd/hoffman2/acme_scratch/edison/20190423.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.edison.restrictedMelt/run']
meshpath = ['/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc',
            '/project/projectdirs/acme/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc']

vartitle = ['z','T','S','rho','u','v','U','taux','tauy','tau','ssh','curl','zice']
varlabel = ['Depth (m)','Temp (degC)','S (PSU)','$\rho$','Velocity +East (m/s)','Velocity +North (m/s)',
            'Velocity magnitude (m/s)','Surface stress','Sea Surface Height (m)','curl']
surfvar = ['taux','tauy','tau','ssh']
varname  = ['depth','timeMonthly_avg_activeTracers_temperature',
            'timeMonthly_avg_activeTracers_salinity',
            'timeMonthly_avg_potentialDensity',
            'timeMonthly_avg_velocityZonal',
            'timeMonthly_avg_velocityMeridional',
            '','timeMonthly_avg_windStressZonal',
            'timeMonthly_avg_windStressMeridional',
            '','timeMonthly_avg_ssh','landIceDraft']

varmin = [str() for i in vartitle]
varmax = [str() for i in vartitle]

# variable axis limits 
# choices for on the continental shelf
varmin = [-1800,-2.5, 34.2, 1027.40, -0.02, -0.02, -0.02,-0.003,-0.003,-0.003,-5.]
varmax = [-100,-0.5, 34.6, 1027.75,  0.02,  0.02, 0.02,0.003,0.003,0.003,-2.]
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
varcmap = ['cmo.deep','cmo.thermal','cmo.haline','cmo.dense','cmo.balance','cmo.balance','cmo.speed','cmo.balance','cmo.balance','cmo.speed','cmo.deep','cmo.curl']
  
#----------------------------------------------------------------------
# TSERIES1 
# -- plot timeseries of variable at a given depth and geographic location
# -- plot geographic location on an area map with bathymetry
#
# Inputs:
#   run       runname, string
#   varlist   variables to plot, list of strings
#   latS      latitude, always in Southern Hem, positive, real 
#   lonW      longitude, always in Western Hem, positive, real
#   startyr   lower limit on simulated year to plot, real
#   endyr     upper limit on simulated year to plot, real 
#   z         depth value, real
#   zab       true if z denotes m above sea level, false if m below surface
#   runcmp    if true, plot both entries in runname
#   savepath  path to save plot images
#----------------------------------------------------------------------
def tseries1(run,varlist,latS,lonW,startyr,endyr,z=0,zab=False,
             runcmp=False,savepath='/global/homes/c/cbegeman/weddell_output/'):

    years = np.arange(startyr,endyr+1,1)
    months = np.arange(1,13,1)
    nt = len(years)*len(months)
    times = np.zeros((nt,))
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    latCell = fmesh.variables['latCell'][:]
    lonCell = fmesh.variables['lonCell'][:]
    xCell = fmesh.variables['xCell'][:]
    yCell = fmesh.variables['yCell'][:]
    locname = str(latS) + 'S' + str(lonW) + 'W'
    latplt = -1.*latS*deg2rad
    lonplt = (360.-lonW)*deg2rad
    idx = np.argmin( (latCell-latplt)**2 + (lonCell-lonplt)**2)
    
    # northern limit for subplots
    logical_N = (latCell < lat_N*deg2rad) & (xCell > 0)
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    zice     = fmesh.variables['landIceDraft'][0,:]
    
    placename = ( str(latS) + 'S' + str(lonW) + 'W' ) 

    if not os.path.exists(savepath + 'bathy_' + placename + '.png'):
        fig = plt.figure()
        cntr1 = plt.tricontourf(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                                zmax[logical_N].flatten(), 20, cmap="viridis")
        plt.plot(yCell[logical_N],     xCell[logical_N], 'o', color = 'white', 
                 markersize = 4, fillstyle = 'none')#, alpha = 0.5)
        plt.plot(yCell[idx],     xCell[idx], 'o', color = 'red', 
                 markersize = 4)
        cntr = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                               zice[logical_N].flatten(), [-10], colors = 'k')
        if varlim:
            plt.clim([varmin[vartitle.index('z')], varmax[vartitle.index('z')]])
        plt.axis('equal')
        cbar = plt.colorbar(cntr1)
        cbar.set_label('Depth (m)')    
        plt.savefig(savepath + 'bathy_' + placename + '.png',dpi=set_dpi)
        plt.close()
    
    data = np.zeros((len(varlist),nt)) 
    if runcmp:
        data2 = np.zeros((len(varlist),nt)) 
    
    t=0
    for yr in years:
        for mo in months:
            times[t] = yr+(mo-1.0)/12.0
       
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
            f = netCDF4.Dataset(filename, 'r')
            for i,var in enumerate(varlist):
                if var not in surfvar:
                    zmid = zmidfrommesh(fmesh,zidx=[idx])
                    if zab:
                        zeval = np.subtract(zmid[-1],z)
                    else:
                        zeval = -1*z
                    zidx = np.argmin(np.abs(np.subtract(zmid,zeval)))
 
                if var in surfvar:
                    data[i,t] = f.variables[varname[vartitle.index(var)]][0,idx]
                else:
                    data[i,t] = f.variables[varname[vartitle.index(var)]][0,idx,zidx]
            f.close()
            if runcmp:
                filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index('ISMF-noEAIS')]) + datestr + '-01.nc'
                f = netCDF4.Dataset(filename, 'r')
                for i,var in enumerate(varlist):
                    if var in surfvar:
                        data2[i,t] = f.variables[varname[vartitle.index(var)]][0,idx]
                    else:
                        data2[i,t] = f.variables[varname[vartitle.index(var)]][0,idx,zidx]
                f.close()    
            t=t+1
    
#    size = 5
#    fig = plt.figure()
#    sc1=plt.scatter(data[0,:],data[1,:],s=size,c=times,cmap='jet')
#    fig.colorbar(sc1)
#    plt.savefig(savepath + run + '_' + varlist[0] + varlist[1] + '_sc_' + placename + '_' + str(startyr) + '-' + str(endyr) + '.png')
#    plt.close()

    nrow=len(varlist)
    ncol=1
    fig,axvar = plt.subplots(nrow,ncol,sharex=True)
    plt.title(run + ': ' + locname)
    for i,var in enumerate(varlist):
        if i == nrow-1:
            axvar[i].set(xlabel='year')
        axvar[i].set(ylabel=var)
        pc = axvar[i].plot(times,data[i,:],'-k')
        if runcmp:
            pc2 = axvar[i].plot(times,data2[i,:],'-b')
    
    if zab:
        m = 'mab'
    else:
        m = 'm'
    filename = run + '_' + varlist[0] + varlist[1] + '_t_' + str(z) + m + '_' + placename + '_' + str(startyr) + '-' + str(endyr) 
    plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)

#----------------------------------------------------------------------
# HOVMOLLER 
# -- color plot of variable vs. depth and time 
#
# Inputs:
#   run       runname, string
#   latS      latitude, always in Southern Hem, positive, real 
#   lonW      longitude, always in Western Hem, positive, real
#   startyr   lower limit on simulated year to plot, real
#   endyr     upper limit on simulated year to plot, real 
#   varlist   variables to plot, list of strings
#   maxDepth  maximum depth of plots
#   savepath  path to save plot images
#----------------------------------------------------------------------
def hovmoller(run,latS,lonW,startyr,endyr,
              varlist = ['T','S','rho','u','v'],maxDepth=-800.0,
              savepath='/global/homes/c/cbegeman/weddell_output/'):

    years = np.arange(startyr,endyr+1,1)
    months = np.arange(1,13,1)
    nt = len(years)*len(months)
    times = np.zeros((nt,))
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    latCell = fmesh.variables['latCell'][:]
    lonCell = fmesh.variables['lonCell'][:]
    xCell = fmesh.variables['xCell'][:]
    yCell = fmesh.variables['yCell'][:]
    locname = str(latS) + 'S' + str(lonW) + 'W'
    latplt = -1.*latS*deg2rad
    lonplt = (360.-lonW)*deg2rad
    idx = np.argmin( (latCell-latplt)**2 + (lonCell-lonplt)**2)
    
    # calculate z from depths
    zmid = zmidfrommesh(fmesh,zidx=[idx])
    z = zmid[0,:]
    nz = fmesh.variables['maxLevelCell'][idx]

    data = np.zeros((len(varlist),nz,nt)) 
    
    t=0
    for yr in years:
        for mo in months:
            times[t] = yr+(mo-1.0)/12.0
       
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
            f = netCDF4.Dataset(filename, 'r')
            for i,var in enumerate(varlist):
                data[i,:,t] = f.variables[varname[vartitle.index(var)]][0,idx,:nz]
            
            f.close()
            t += 1
    
    nrow=len(varlist)
    ncol=1
    fig,axvar = plt.subplots(nrow,ncol,sharex=True)
    #plt.title(run + ': ' + locname)
    
    for i,var in enumerate(varlist):
    
        if i == nrow-1:
            axvar[i].set(xlabel='year')
        axvar[i].set(ylabel=var)
        pc = axvar[i].pcolor(times, z[:nz], data[i,:], 
                        vmin=varmin[vartitle.index(var)], vmax=varmax[vartitle.index(var)], 
                        cmap=varcmap[vartitle.index(var)])
        axvar[i].set_ylim((maxDepth, 0))
        fig.colorbar(pc,ax=axvar[i])
    
    filname = run + '_hovmoller_' + varlist[0] + varlist[1] + locname + '_' + str(startyr) + '-' + str(endyr)
    plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
    plt.close()

#----------------------------------------------------------------------
# PROFILE
# -- variable vs. depth at a given geographic location
#
# Inputs:
#   varlist   variables to plot, list of strings
#   run       runname, string
#   startyr   lower limit on simulated year to plot, real
#   endyr     upper limit on simulated year to plot, real 
#   latS      latitude, always in Southern Hem, positive, real 
#   lonW      longitude, always in Western Hem, positive, real
#   maxDepth  maximum depth of plots
#   runcmp    if true, plot both entries in runname
#   mo        month to plot data from, if 0 then plot all months of data
#   savepath  path to save plot images
#----------------------------------------------------------------------
def profile(varlist,run,startyr,endyr,latS,lonW,
            maxDepth = -500.,runcmp=False,mo = 0,
            savepath='/global/homes/c/cbegeman/weddell_output/'):

    varmin[vartitle.index('T')] = -2.2
    varmax[vartitle.index('T')] = 2
    varmin[vartitle.index('S')] = 32.8
    varmax[vartitle.index('S')] = 34.8
    varmin[vartitle.index('rho')] = 1026.7
    varmax[vartitle.index('rho')] = 1027.8
    varmin[vartitle.index('u')] = -0.03
    varmax[vartitle.index('u')] = 0.03
    varmin[vartitle.index('v')] = -0.03
    varmax[vartitle.index('v')] = 0.03
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    
    latCell = fmesh.variables['latCell'][:]
    lonCell = fmesh.variables['lonCell'][:]
    xCell    = fmesh.variables['xCell'][:]
    yCell    = fmesh.variables['yCell'][:]
    depths = fmesh.variables['refBottomDepth'][:]
    z = np.zeros(depths.shape)
    z[0] = -0.5 * depths[0]
    z[1:] = -0.5 * (depths[0:-1] + depths[1:])
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    zice     = fmesh.variables['landIceDraft'][0,:]
    logical_N = (latCell < lat_N*deg2rad) & (xCell > 0)
    
    locname = str(latS) + 'S' + str(lonW) + 'W'
    latplt = -1.*latS*deg2rad
    lonplt = (360.-lonW)*deg2rad
    idx = np.argmin( (latCell-latplt)**2 + (lonCell-lonplt)**2)  #122901-1
    
    if not os.path.exists(savepath + 'bathy_' + locname + '.png'):
        fig = plt.figure()
        cntr1 = plt.tricontourf(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                                zmax[logical_N].flatten(), 20, cmap="viridis")
        plt.plot(yCell[logical_N],     xCell[logical_N], 'o', color = 'white', 
                 markersize = 4, fillstyle = 'none')#, alpha = 0.5)
        plt.plot(yCell[idx],     xCell[idx], 'o', color = 'red', 
                 markersize = 4)
        cntr = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                               zice[logical_N].flatten(), [-10], colors = 'k')
        plt.axis('equal')
        cbar = plt.colorbar(cntr1)
        cbar.set_label('Depth (m)')    
        plt.savefig(savepath + 'bathy_' + locname + '.png',dpi=set_dpi)
        plt.close()

    years = np.arange(startyr,endyr+1,1)
    if mo == 0:
       months = np.arange(1,13,1)
    else:
       months = [mo]
    nt = len(years)*len(months)
    times = np.zeros((nt,))
    colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, 13)]
    
    lineStyle = ['-',':','-.','--']
    
    nrow=1
    ncol=len(varlist)
    fig,axvar = plt.subplots(nrow,ncol,sharey=True)
    axvar[0].set(ylabel='depth (m)')
    
    t = 0
    for i,yr in enumerate(years):
        for j,mo in enumerate(months):
            c = colors[j]
            
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
            f = netCDF4.Dataset(filename, 'r')
            for k,var in enumerate(varlist):
            
               #if i > 0:
               #    axvar[i] = fig.add_subplot(nrow, ncol, 2, sharey=axvar[0])
               #else:
               #    axvar[i] = fig.add_subplot(nrow, ncol, i+1)
               axvar[k].set(xlabel=var)
            
               data     = f.variables[varname[vartitle.index(var)]]
               if runcmp:
                  filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index('ISMF-noEAIS')]) + datestr + '-01.nc'
                  f2 = netCDF4.Dataset(filename, 'r')
                  data2     = f2.variables[varname[vartitle.index(var)]]
                  f2.close()
               f.close()    
               axvar[k].plot(data[0,idx,:], z, label="yr{0:04d}".format(yr), 
                             color=c,linestyle=lineStyle[i])
               if runcmp:
                  axvar[k].plot2(data2[0,idx,:], z, label="yr{0:04d}".format(yr), 
                                color=c,linestyle='--')
    
               axvar[k].set_xlim([varmin[vartitle.index(var)], varmax[vartitle.index(var)]])
               axvar[k].set_ylim([maxDepth,0])
               axvar[k].grid()
    datestr=str(startyr) + '-' + str(endyr)
    plt.title(run + ': ' + datestr)
    filename = savepath + run 
    if runcmp:
       filename = filename + '_cmp'
    filename = filename + '_profiles_' + varlist[0] + varlist[1] + '_' + locname + '_' + datestr
    print(filename)
    plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
    plt.close()

#----------------------------------------------------------------------
# TRANSECT
# -- plot variable in color vs. depth across a profile in lon or lat
#
# Inputs:
#   latS      latitude range, always in Southern Hem, positive, 
#             vector of length 2, real 
#   lonW      longitude range, always in Western Hem, positive,
#             vector of length 2, real 
#   varlist   variables to plot, list of strings
#   yr   lower limit on simulated year to plot, real
#   mo        month to plot data from, if 0 then plot all months of data

# optional variables:
#   zscale    scale of depth axis
#   run       runname, string 
#   runcmp    if true, plot difference between variables of both entries in runname
#   savepath  path to save plot images
#----------------------------------------------------------------------
def transect(latS, lonW, varlist, yr, mo, varlim = False, zscale = 'linear', 
             run='ISMF', runcmp = False, new = False,
             savepath='/global/homes/c/cbegeman/weddell_output/'):

    placename = ( str(latS[0]) + 'S' + str(lonW[0]) + 'W-' + 
                  str(latS[1]) + 'S' + str(lonW[1]) + 'W' )
    
    lat_transect = np.multiply(-1.,latS)
    lon_transect = np.subtract(360,lonW)
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    datestr = '{0:04d}-{1:02d}'.format(yr, mo)
    filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
    f = netCDF4.Dataset(filename, 'r')
    if runcmp:
        filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index('ISMF-noEAIS')]) + datestr + '-01.nc'
        f2 = netCDF4.Dataset(filename, 'r')

    # constants
    dlat = 0.15 # at 30km resolution, distance between cells in latitude space
    dlon = 0.98
    
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
    for var in varlist:
        if os.path.exists(savepath + run + '_' + var + '_' + placename + '_' + datestr + '_lim' + str(varlim)+'.png'):
            print(run + '_' + var + '_' + placename + '_' + datestr + '_lim' + str(varlim) + ' exists')
            #continue
    
        data     = f.variables[varname[vartitle.index(var)]][0,:]
        if runcmp: 
            data2 = f2.variables[varname[vartitle.index(var)]][0,:]
        
        # calculate z from depths
        zmid = zmidfrommesh(fmesh)
        #zssh = zbottom[:,0] + zh[:,0]
        
        # define plot location
        if new:
            idx1 = np.argmin( np.square(yCell-lat_transect[0]) + 
                              np.square(xCell-lon_transect[0])   )
            idx2 = np.argmin( np.square(yCell-lat_transect[1]) + 
                              np.square(xCell-lon_transect[1])   )
            return
        else:
            # define line of constant latitude
            if lat_transect[0] == lat_transect[1]:
                lat_transect[0] = lat_transect[0] - dlat
                lat_transect[1] = lat_transect[1] + dlat
            if lon_transect[0] == lon_transect[1]:
                lon_transect[0] = lon_transect[0] - dlon
                lon_transect[1] = lon_transect[1] + dlon
            # northern limit for subplots
            logical_N = (latCell < lat_N*deg2rad) & (xCell > 0)
            
            # indices of transect
            logical_trans = ( (latCell > lat_transect[0]*deg2rad) & 
                              (latCell < lat_transect[1]*deg2rad) &
                              (lonCell > lon_transect[0]*deg2rad) & 
                              (lonCell < lon_transect[1]*deg2rad)   )
            idx_trans = np.where(logical_trans)[0]
            
            idx1 = np.argmin(yCell[logical_trans])
            temp = np.sqrt( np.square(yCell[logical_trans] - yCell[logical_trans][idx1]) + 
                            np.square(xCell[logical_trans] - xCell[logical_trans][idx1])   )
            idxsort_trans = idx_trans[temp.argsort()]
            ysort_trans = yCell[idxsort_trans]
            xsort_trans = xCell[idxsort_trans]
            dist_trans = temp[temp.argsort()]
            dd_trans = np.zeros(dist_trans.shape)
            dd_trans[1:] = dist_trans[1:]-dist_trans[:-1]
            idx2, = np.where(dd_trans>100e3) 
            if len(idx2) > 0:
                dmax = dist_trans[idx2[0]-1]# distance along transect
            else:
                dmax = np.max(dist_trans)
       
        # create mesh variables for plotting
        ymesh,temp = np.meshgrid(dist_trans, depths);
        zmesh = np.transpose(zmid[idxsort_trans])
        data_trans = np.transpose(data[idxsort_trans,:])
        data_trans_masked = np.transpose(data[idxsort_trans,:])
        data_trans_zmasked = np.transpose(data[idxsort_trans,:])
        # Ideally, mask both bad data and ice shelf
        #data_trans_masked = np.ma.masked_where( (data_trans < bad_data) |
        #                                        (zmesh > zicemesh), data_trans)
        data_trans_masked = np.ma.masked_where( (data_trans < bad_data) |
                (data_trans == bad_data2), data_trans)
        for idx,i in enumerate(idxsort_trans):
            data_trans_zmasked[:,idx] = np.ma.masked_where( 
                    (zmid[i,:] < zmid[i,kmax[i]-1]), data[i,:])
        mask = np.ma.getmask(data_trans_masked)
        
        temp1 = kmax[logical_trans]
        temp2 = kmax[idxsort_trans]
        # plots
        
        # show profile line across cells
        if not os.path.exists(savepath + 'bathy_' + placename + '.png'):
            fig1 = plt.figure()
            plt.plot(yCell[logical_N],     xCell[logical_N], 'k.')
            plt.plot(yCell[logical_trans], xCell[logical_trans], 'r.')
            plt.axis('equal')
            plt.savefig('grid_' + placename + '_' + datestr + '.png',dpi=set_dpi)
            plt.close()
            
            fig = plt.figure()
            cntr1 = plt.tricontourf(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                                    zmax[logical_N].flatten(), 20, cmap="viridis")
            plt.plot(yCell[logical_N],     xCell[logical_N], 'o', color = 'white', 
                     markersize = 4, fillstyle = 'none')#, alpha = 0.5)
            cntr = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                                   zice[logical_N].flatten(), [-10], colors = 'k')
            plt.plot(yCell[idxsort_trans], xCell[idxsort_trans], 'k-')
            plt.axis('equal')
            cbar = plt.colorbar(cntr1)
            cbar.set_label('Depth (m)')    
            plt.savefig(savepath + 'bathy_' + placename + '.png',dpi=set_dpi)
            plt.close()
        
        yline = np.divide(dist_trans,1e3)
        yfill = np.append(yline[0],yline)
        yfill = np.append(yfill,yline[-1])
        sshfill = np.append(0,ssh[idxsort_trans])
        sshfill = np.append(sshfill,0)
        bathymax = np.min(zmax[idxsort_trans]) - 100
        bathyfill = np.append(bathymax,zmax[idxsort_trans])
        bathyfill = np.append(bathyfill,bathymax)
        
        
        clevels = np.arange(varmin[vartitle.index(var)], 
                           varmax[vartitle.index(var)],
                           dvar[vartitle.index(var)]   )
        if clevels[0] > np.min(data_trans_zmasked[~mask].flatten()):
            clevels = np.append(np.min(data_trans_zmasked[~mask].flatten()),clevels)
        if clevels[-1] < np.max(data_trans_zmasked[~mask].flatten()):
            clevels = np.append(clevels,np.max(data_trans_zmasked[~mask].flatten()))
        
        fig2 = plt.figure()
        if runcmp:
            cmap1 = "cmo.balance"
        else:
            cmap1 = varcmap[vartitle.index(var)]
        if runcmp:
            data1_trans_zmasked = data_trans_zmasked
            data2_trans_zmasked = np.transpose(data2[idxsort_trans,:])
            data_trans_zmasked = np.subtract(data1_trans_zmasked,data2_trans_zmasked)
            cntr2 = plt.tricontourf(np.divide(ymesh[~mask].flatten(),1e3), np.abs(zmesh[~mask].flatten()), 
                                data_trans_zmasked[~mask].flatten(), 
                                levels=clevels,cmap=cmap1)
        else:
            cntr2 = plt.tricontourf(np.divide(ymesh[~mask].flatten(),1e3), np.abs(zmesh[~mask].flatten()), 
                                data_trans_zmasked[~mask].flatten(), 
                                levels=clevels,cmap=cmap1)
        plt.plot(yline, np.abs(ssh[idxsort_trans]), color = 'black', marker = '.', linestyle = '-')
        plt.plot(np.divide(ymesh[~mask].flatten(),1e3), np.abs(zmesh[~mask].flatten()), 
                 '.', color = 'white', markersize = 1)#, fillstyle = 'none')
        plt.fill(yfill, np.abs(sshfill), c = 'white', alpha = 1)
        plt.plot(np.divide(dist_trans,1e3), np.abs(zmax[idxsort_trans]), color = 'black', marker = '.', linestyle = '-')
        plt.fill(yfill, np.abs(bathyfill), c = 'grey', alpha = 1)
        ax = plt.gca()
        if varlim:
            if runcmp:
                if var == 'T':
                    cmax = 1.
                else:
                    cmax = max(abs(np.percentile(data_trans_zmasked[~mask].flatten(),10)),
                           abs(np.percentile(data_trans_zmasked[~mask].flatten(),90)))
                plt.clim([-1*cmax,cmax])
            else:
                plt.clim([varmin[vartitle.index(var)], varmax[vartitle.index(var)]])
        if zscale == 'log':
            ax.set_yscale('log')
            plt.ylim([1,np.abs(np.min(zmesh[~mask].flatten())-50)])
        ax.invert_yaxis()
        cbar = plt.colorbar()
        cbar.set_label(var)
        plt.xlabel('Distance (km)')
        plt.ylabel('Depth (m)')
        plt.xlim([np.min(dist_trans)/1e3,dmax/1e3])
        plt.title(run + ': ' + datestr)
        filename = run
        if runcmp:
            filename = filename + '_cmp'
        filename = filename + '_' + var + '_' + placename + '_' + datestr + '_lim' + str(varlim) 
        plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
        print(filename) 
        plt.close()

#------------------------------------------------------------------------------
# PLOT_MESH_VAR 
# -- Plots a quantity in map view and saves figure to file
#
# Inputs:
# optional variables:
#    run        name of the model run, string
#    locname    name of the location which sets map limits, string
#    savepath   directory where figure is to be saved, string
#------------------------------------------------------------------------------
def plot_zice_map(run = 'ISMF',locname = 'fris',savepath='/global/homes/c/cbegeman/weddell_output/'):

    # set further parameters
    lat_N = -50 # northern limit of domain
    
    if locname not in loc:
        print('locname is not defined')
        return
    xmax = loc_xmax[loc.index(locname)]
    ymin = loc_ymax[loc.index(locname)]
    size = loc_ptsize[loc.index(locname)]
    
    filename = 'fmesh_zice_' + locname
    print(filename)
    if os.path.exists(savepath + filename + '.png'):
        print('file exists')
        if not overwrite:
            return
    
    # open data files
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    
    # import variables from file
    latCell  = fmesh.variables['latCell'][:]
    lonCell  = fmesh.variables['lonCell'][:]
    idxCell  = fmesh.variables['indexToCellID'][:]
    xCell    = fmesh.variables['xCell'][:]
    yCell    = fmesh.variables['yCell'][:]
    kmax     = fmesh.variables['maxLevelCell'][:]
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    zice     = fmesh.variables['landIceDraft'][0,:]
    
    # northern limit for subplots
    idx_N = np.argwhere((xCell<xmax) & (xCell > 0) & (yCell < 0) & (yCell > -1.5e6) & (latCell < lat_N*deg2rad))
    logical_N = idx_N[:,0]
    x = xCell[logical_N]
    y = yCell[logical_N]
    idx = idxCell[logical_N]

    fig = plt.figure()
    
    ax = fig.gca()
    ax.set_aspect('equal')
    gl2 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                         zmax[logical_N].flatten(), [-1800], colors = 'k', linewidths = 2)
    gl1 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                         zice[logical_N].flatten(), [-10], colors = 'b', linewidths = 2)
    
    cntr1 = plt.scatter(y, x, s=size, c=zice[logical_N], marker = '.',
                        cmap = "cmo.deep") 
    
    cbar1 = fig.colorbar(cntr1)
    #cbar1.set_label(var)
    
    plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
    return

#------------------------------------------------------------------------------
# plot_surf_var
#
# Description: Plots a quantity in map view and saves figure to file
#
# Inputs:
#    var        variable to be plotted
#    yr         plot at year yr
#    mo         plot at month mo
#    run        (optional) name of run. Default ISMF
#    locname    location name (string) sets geographic bounds on plot
#    plottype   option for velocity fields (string):
#               'abs' magnitude of velocity vector (color plot)
#               'quiver' quiver plot of velocity field
#               'head' orientation of velocity vector (color plot)
#    z          (optional) meters below surface to plot var. z is negative 
#               below the surface and 0 at the surface. Default is 0.
#    zab        (optional) indicates that z is defined as meters above bottom 
#               (z must be >0)
#    level      (optional) Indicates the value of var at which to plot contour 
#               depth
#    runcmp     (optional) logical value to indicate whether two runs should be 
#               differenced
#    varlim     (optional) logical value to set color limits for var
#    overwrite  if true, create plot even if it already exists
#               if false, skip all [var,yr,mo] combinations for which plot already exists
#    savepath   directory where figure is to be saved
#------------------------------------------------------------------------------
def plot_surf_var(var,yr,mo,run='ISMF',locname='fris',plottype = 'abs',
                  z=0,zab=False,level=bad_data,
                  runcmp=False,varlim=False,overwrite=False,
                  savepath='/global/homes/c/cbegeman/weddell_output/'):

    if locname not in loc:
        print('locname is not defined')
        return
    xmax = loc_xmax[loc.index(locname)]
    ymin = loc_ymax[loc.index(locname)]
    size = loc_ptsize[loc.index(locname)]
    
    if zab:
        m = 'mab'
    else:
        m = 'm'
    
    # variable is only defined at the surface
    if (var in surfvar):
        z = 0.
    
    datestr = '{0:04d}-{1:02d}'.format(yr, mo)
    
    # filename to save plot
    filename = run
    if runcmp:
        filename = filename + '_cmp'
    filename = filename + '_' + var + '_' 
    if plottype == 'quiver':
        filename = filename + 'quiver_'
    elif plottype == 'head':
        filename = filename + 'heading_'
    filename = filename + locname + '_'
    if level == bad_data:
        filename = filename + str(int(abs(z))) + m 
    else:
        print(level)
        if var == 'rho':
            filename = filename + str(round(abs((level-1000)*10))) + 'c'
        else:
            filename = filename + str(round(abs(level))) + 'c'
    filename = filename + '_' + datestr 
    if varlim:
        filename = filename + '_limTrue'
    print(filename)
    
    # check if plot was already generated
    if os.path.exists(savepath + filename + '.png'):
        print('file exists')
        if not overwrite:
            print('skipping file')
            return
    
    # open data files
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    filein = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
    f = netCDF4.Dataset(filein, 'r')
    if runcmp:
       filein2 = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index('ISMF-noEAIS')]) + datestr + '-01.nc'
       f2 = netCDF4.Dataset(filein2, 'r')
    
    # import variables from file
    latCell  = fmesh.variables['latCell'][:]
    lonCell  = fmesh.variables['lonCell'][:]
    idxCell  = fmesh.variables['indexToCellID'][:]
    xCell    = fmesh.variables['xCell'][:]
    yCell    = fmesh.variables['yCell'][:]
    kmax     = fmesh.variables['maxLevelCell'][:]
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    zh       = fmesh.variables['layerThickness'][0,:]
    zice     = fmesh.variables['landIceDraft'][0,:]
    
    # northern limit for subplots
    idx_N = np.argwhere((xCell<xmax) & (xCell > 0) & (yCell < 0) & (yCell > -1.5e6) & (latCell < lat_N*deg2rad))
    logical_N = idx_N[:,0]
    x = xCell[logical_N]
    y = yCell[logical_N]
    idx = idxCell[logical_N]
    
    # calculate z from depths
    zmid = zmidfrommesh(fmesh)
    
    # get data at specified depth
    zidx = np.zeros(len(logical_N),)
    # define the index for specifed depth
    if level == bad_data:
        if zab:
            if z == 0:
                zidx = kmax[logical_N]-1
            else:
                zeval = np.add(zmax[logical_N],z)
                for i,idx in enumerate(logical_N):
                    zidx[i] = int(np.argmin(np.abs(np.subtract(zmid[idx,:],zeval))))
        else:
            for i,idx in enumerate(logical_N):
                zidx[i] = int(np.argmin(np.abs(np.subtract(zmid[idx,:],-1*z))))
        
        # get data
        if var == 'U':
            u = f.variables['timeMonthly_avg_velocityZonal'][0,logical_N,:]
            v = f.variables['timeMonthly_avg_velocityMeridional'][0,logical_N,:]
            uz = np.zeros(len(zidx),)
            vz = np.zeros(len(zidx),)
            for i,idx in enumerate(zidx):
                uz[i] = u[i,int(idx)]
                vz[i] = v[i,int(idx)]
            bad_idx = (uz > bad_data) #| (heading != bad_data2)
            dataz = np.sqrt(np.add(np.square(uz),np.square(vz)))
            heading = np.divide(np.arctan2(uz,vz),deg2rad)
        elif var == 'tau':
            uz = f.variables['timeMonthly_avg_windStressZonal'][0,logical_N]
            vz = f.variables['timeMonthly_avg_windStressMeridional'][0,logical_N]
            bad_idx = (uz > bad_data) #| (heading != bad_data2)
            dataz = np.sqrt(np.add(np.square(uz),np.square(vz)))
            heading = np.divide(np.arctan2(uz,vz),deg2rad)
        elif var == 'ssh':
            dataz = f.variables[varname[vartitle.index(var)]][0,logical_N]
            bad_idx = (dataz > bad_data) #| (heading != bad_data2)
        else:
            data = f.variables[varname[vartitle.index(var)]][0,logical_N,:]
            if (var not in surfvar):
                dataz = np.zeros(len(zidx),)
                for i,idx in enumerate(zidx):
                    dataz[i] = data[i,int(idx)]
            else:
                dataz = data
            bad_idx = (dataz > bad_data) #| (heading != bad_data2)
    
    # get the depth at which data equals value specficied in level
    else:
        data = f.variables[varname[vartitle.index(var)]][0,logical_N,:]
        dataz = np.zeros(len(zidx),)
        icemask = np.zeros(len(zidx),dtype=bool)
        sfmask = np.zeros(len(zidx),dtype=bool)
        for i in range(0,len(zidx),1):
            zidx = np.argmin(np.abs(np.subtract(data[i,:],level)))
            cells,zlevels = zh.shape
            dataz[i] = zmid[i,zidx]
            if (zidx == 0 or min(data[i,:])>level):
                icemask[i] = True
            if (zidx == zlevels-1 or max(data[i,:])<level):
                sfmask[i] = True
            #dataz[i] = max(zmid[i,:])
        bad_idx = (dataz > bad_data) #| (heading != bad_data2)
        dataz = dataz + zice[logical_N]
        dataz = np.multiply(dataz,-1.)
 
    # solve for difference between fields across model runs
    if runcmp:
        dataz1 = dataz
        if var == 'U':
            heading1 = heading
            uz1 = uz
            vz1 = vz
            u = f2.variables['timeMonthly_avg_velocityZonal'][0,logical_N,:]
            v = f2.variables['timeMonthly_avg_velocityMeridional'][0,logical_N,:]
            uz2 = np.zeros(len(x),)
            vz2 = np.zeros(len(x),)
            i = 10
            for i,idx in enumerate(zidx):
                uz2[i] = u[i,int(idx)]
                vz2[i] = v[i,int(idx)]
            #for i in range(0,len(x),1):
            #    uz2[i] = u[i,zidx]
            #    vz2[i] = v[i,zidx]
            dataz2 = np.sqrt(np.add(np.square(uz2),np.square(vz2)))
            heading2 = np.divide(np.arctan2(uz2,vz2),deg2rad)
            dataz = np.subtract(dataz1,dataz2)
            heading = np.subtract(heading1,heading2)
            uz = np.subtract(uz1,uz2)
            vz = np.subtract(vz1,vz2)
        elif var == 'tau':
            heading1 = heading
            uz1 = uz
            vz1 = vz
            uz2 = f2.variables['timeMonthly_avg_windStressZonal'][0,logical_N]
            vz2 = f2.variables['timeMonthly_avg_windStressMeridional'][0,logical_N]
            bad_idx = (uz2 > bad_data) #| (heading != bad_data2)
            dataz2 = np.sqrt(np.add(np.square(uz2),np.square(vz2)))
            heading2 = np.divide(np.arctan2(uz2,vz2),deg2rad)
            dataz = np.subtract(dataz1,dataz2)
            heading = np.subtract(heading1,heading2)
            uz = np.subtract(uz1,uz2)
            vz = np.subtract(vz1,vz2)
        elif var == 'ssh':
            dataz2 = f2.variables['timeMonthly_avg_ssh'][0,logical_N]
            dataz = np.subtract(dataz1,dataz2)
        else:
            data2 = f2.variables[varname[vartitle.index(var)]][0,logical_N,:]
            if (var not in surfvar):
                dataz2 = np.zeros(len(data2),)
                print('load data2 for clevel')
                print(data2.shape)
                for i in range(0,len(data2),1):
                    if level > bad_data:
                        zidx = np.argmin(np.abs(np.subtract(data2[i,:],level)))
                        dataz2[i] = zmid[i,int(zidx)]
                    else:
                        dataz2[i] = data2[i,int(zidx[i])]
            else:
                dataz2 = data2
            bad_idx = (dataz2 > bad_data) #| (heading != bad_data2)
            dataz = np.subtract(dataz1,dataz2)
    
    # plots
    if plottype == 'abs':
        fig = plt.figure()
        
        if runcmp:
            cmap1 = "cmo.balance"
        else:
            cmap1 = varcmap[vartitle.index(var)]
            if level != bad_data:
                 cmap1 = "cmo.deep"
            #if runcmp:
            #    cntr1 = plt.scatter(y[bad_idx], x[bad_idx], 
            #                s = size, c = dataz[bad_idx], marker = '.',
            #                cmap = cmap1)#,vmin=-0.01,vmax=0.01)
            #else:
        cntr1 = plt.scatter(y[bad_idx], x[bad_idx], 
                    s = size, c = dataz[bad_idx], marker = '.',
                    cmap = cmap1) #, norm = LogNorm())
        if varlim:
            plt.clim([varmin[vartitle.index(var)], varmax[vartitle.index(var)]])
        
        elif level > bad_data and not runcmp:
            #plt.clim([0,2000])
            plt.clim([np.percentile(dataz[bad_idx],10),np.percentile(dataz[bad_idx],90)])
        
        if (var == 'U'):
        #if runcmp:
            cmax = max(abs(np.percentile(dataz[bad_idx],10)),abs(np.percentile(dataz[bad_idx],90)))
            plt.clim([-1*cmax,cmax])

        ax = fig.gca()
        ax.set_aspect('equal')
        
        # plot an arbitrarily deep contour as a rough marker of continental shelf break
        gl2 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                             zmax[logical_N].flatten(), [-1800], colors = 'k', linewidths = 2)
        # plot a fairly shallow ice contour as a marker of ice shelf front
        gl1 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                             zice[logical_N].flatten(), [-10], colors = 'b', linewidths = 2)
        
        if level > bad_data and not runcmp:
            # outline points where the contour is below seafloor
            cntro = plt.scatter(y[sfmask], x[sfmask], 
                            s=size, marker = '.',
                            facecolor = 'none',edgecolors = 'k')
            # outline points where the contour is above ice base
            cntro = plt.scatter(y[icemask], x[icemask], 
                            s=size, marker = '.',
                            facecolor = 'none',edgecolors = 'b')
             
        cbar1 = fig.colorbar(cntr1)
        if level > bad_data:
            cbar1.set_label('Contour depth (m)')
        else:
            cbar1.set_label(var)
        fig.tight_layout()
        plt.title(run + ': ' + datestr)
        plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
        plt.close()
    
    if (var == 'U') | (var == 'tau'):
        
        if plottype == 'head':
            fig = plt.figure()
            if runcmp:
                cntr2 = plt.scatter(y[bad_idx], x[bad_idx],
                                s=size, c=heading[bad_idx],  marker = '.', cmap="cmo.balance",vmin=-10,vmax=10) 
            else:
                cntr2 = plt.scatter(y[bad_idx], x[bad_idx],
                                s=size, c=heading[bad_idx],  marker = '.', cmap="cmo.phase", 
                                vmin = -180, vmax = 180)
            gl2 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                                   zice[logical_N].flatten(), [-10], colors = 'k', linewidths = 2)
            ax = fig.gca()
            ax.set_aspect('equal')
            cbar2 = fig.colorbar(cntr2)
            cbar2.set_label('Degrees from N')    
            fig.tight_layout()
            if runcmp:
                run = run + '_cmp'
            plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
            plt.close()
        
        elif plottype == 'quiver':
            fig = plt.figure()
            cntr1 = plt.tricontourf(y.flatten(), x.flatten(), 
                                    zmax[logical_N].flatten(), 20, 
                                    cmap="cmo.deep")
            qvr = plt.quiver(y[bad_idx], x[bad_idx], uz[bad_idx], vz[bad_idx],color='m')
            gl2 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
                                   zice[logical_N].flatten(), [-10], colors = 'b', linewidths = 2)
            if varlim:
                plt.clim([varmin[vartitle.index('z')], varmax[vartitle.index('z')]])
            ax = fig.gca()
            ax.set_aspect('equal')
            cbar2 = fig.colorbar(cntr1)
            cbar2.set_label('depth (m)')    
            fig.tight_layout()
            plt.title(run + ': ' + datestr)
            plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
            plt.close()

#------------------------------------------------------------------------------
# ZMIDFROMMESH 
# -- solve for z values at the middle of the cell where prognostic variables are solved 
#
# Inputs:
#    fmesh  data structure for the mesh
#    zidx   (optional) indices for z level (list of integers)
#           if zidx=-1 (default), return the full z vector
#------------------------------------------------------------------------------
def zmidfrommesh(fmesh,zidx=[-1]):
    zh       = fmesh.variables['layerThickness'][0,:]
    cells,zlevels = zh.shape
    # calculate z from depths
    if zidx[0]<0:
        zidx = np.ones((cells,),dtype=bool)
    zmax     = fmesh.variables['bottomDepth'][:]
    zloc = zmax[zidx]
    zbottom  = np.zeros(zh[zidx,:].shape)
    zbottom[:,-1] = np.multiply(-1,zmax[zidx])
    for i in range(zlevels-2,-1,-1):
        zbottom[:,i] = zbottom[:,i+1] + zh[zidx,i+1]
    zmid = zbottom + np.multiply(0.5,zh[zidx])
    return zmid
