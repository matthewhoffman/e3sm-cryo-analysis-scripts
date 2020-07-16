#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

import sys
sys.path.append('/global/homes/c/cbegeman/e3sm-cryo-analysis-scripts/great_circle_calculator')
import os
import csv
import netCDF4
import cartopy
import pyproj
#import LatLon
import great_circle_calculator.great_circle_calculator as great_circle
#import datetime
import numpy as np
import matplotlib as pltlib
from matplotlib import ticker,rc
import matplotlib.pyplot as plt
import scipy.signal
import scipy.interpolate as interp
import matplotlib.cm as cmx
import matplotlib.colors as colors
from matplotlib.colors import LogNorm,Normalize
from matplotlib.colors import SymLogNorm
import cmocean
import pandas
from extract_depths import zmidfrommesh
from pick_from_mesh import *
from plot_config import *

global bad_data, bad_data2, deg2rad, lat_N, runname, runpath, meshpath, vartitle, varname, varmin, varmax, varcmap, surfvar, dvar

def TS_diagram(run,lat,lon,startyr,endyr,z=0,zab=False,zall=True,plot_lines=True,
               seasonal=False,runcmp=False,savepath=savepath_nersc,
               polygon = TSpolygon_Hattermann2018):

    placename = (str(int(abs(lat))) + 'S' + 
                 str(int(abs(lon-360))) + 'W') 
    
    years = np.arange(startyr,endyr+1,1)
    months = np.arange(1,13,1)
    nt = len(years)*len(months)
    times = np.zeros((nt,))
    
    idx = pick_point(run=run,lat=lat,lon=lon)
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    #latCell = fmesh.variables['latCell'][:]
    #lonCell = fmesh.variables['lonCell'][:]
    #xCell = fmesh.variables['xCell'][:]
    #yCell = fmesh.variables['yCell'][:]
    #latplt = -1.*latS*deg2rad
    #lonplt = (360.-lonW)*deg2rad
    #idx = np.argmin( (latCell-latplt)**2 + (lonCell-lonplt)**2)
    # northern limit for subplots
    #logical_N = (latCell < lat_N*deg2rad) & (xCell > 0)
    #zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    #zice     = fmesh.variables['landIceDraft'][0,:]
    nz = len(fmesh.variables['layerThickness'][0,idx,:]) 
    
    T = np.zeros((nt,nz)) 
    S = np.zeros((nt,nz)) 
    if runcmp:
        T2 = np.zeros((nt,nz)) 
        S2 = np.zeros((nt,nz)) 
    
    zmid,_,_ = zmidfrommesh(fmesh,cellidx=[idx])
    if zab:
        zeval = np.add(zmid[0][-1],z)
        m = 'mab'
    else:
        zeval = -1*z
        m = 'm'
    zidx = np.argmin(np.abs(np.subtract(zmid,zeval)))
   
    t=0
    for yr in years:
        for mo in months:
            times[t] = yr+(mo-1.0)/12.0
       
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
            f = netCDF4.Dataset(filename, 'r')
            T[t,:] = f.variables[varname[vartitle.index('T')]][0,idx,:]
            S[t,:] = f.variables[varname[vartitle.index('S')]][0,idx,:]
            f.close()
            if runcmp:
                filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index('ISMF-noEAIS')]) + datestr + '-01.nc'
                f = netCDF4.Dataset(filename, 'r')
                T2[t,:] = f.variables[varname[vartitle.index(var)]][0,idx,:]
                S2[t,:] = f.variables[varname[vartitle.index(var)]][0,idx,:]
                f.close()    
            t=t+1
   
    fig = plt.figure(1, facecolor='w')
    axTS = fig.add_subplot()
    if zall:
        sc=axTS.scatter(S, T, s=1, c='grey')
    if plot_lines:
        if seasonal:
            cNorm  = Normalize(vmin=0, vmax=1)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='twilight')
            cbtitle = r'Time of Year'
        else:
            cNorm  = Normalize(vmin=startyr, vmax=endyr + 1)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='cmo.deep')
            cbtitle = r'Simulated Year'
        for i,ti in enumerate(times):
            if i > 0:
               if seasonal:
                   colorVal = scalarMap.to_rgba(np.subtract(ti,np.floor(ti)))
               else:
                   colorVal = scalarMap.to_rgba(ti)
               #scz=axTS.plot([S[i-1,zidx],S[i,zidx]], [T[i-1,zidx],T[i,zidx]], 
               #              '-', color=colorVal,linewidth=1)
               Sline = [S[i-1,zidx],S[i,zidx]]
               Tline = [T[i-1,zidx],T[i,zidx]]
               scz=axTS.arrow(Sline[0], Tline[0], Sline[1]-Sline[0], Tline[1]-Tline[0],
                              color=colorVal,linewidth=1)

    else:
        scz=axTS.plot(S[:,zidx], T[:,zidx], s = 30,edgecolor='k', cmap='twilight')
    S_polygon = [polygon[i][0] for i in np.arange(0,4)]
    T_polygon = [polygon[i][1] for i in np.arange(0,4)]
    S_polygon = np.append(S_polygon,polygon[0][0])
    T_polygon = np.append(T_polygon,polygon[0][1])
    print(S_polygon)
    print(T_polygon)
    for i in np.arange(0,4):
        axTS.plot([S_polygon[i],S_polygon[i+1]],
                  [T_polygon[i],T_polygon[i+1]],
                  color='k')
    #axTS.fill(S_polygon, T_polygon,closed=True,
    #          fill = False,edgecolor='k')

    plt.colorbar(scalarMap,label=cbtitle)
    lncolor = 'black'
    lw1 = 1 
    plt.plot([34.0, 34.0], [-1.85, 0.2],  ':', color=lncolor, linewidth=lw1)
    plt.plot([34.5, 34.5], [-1.86, -1.5], ':', color=lncolor, linewidth=lw1)
    plt.plot([34.0, 35.0], [-1.5, -1.5],  ':', color=lncolor, linewidth=lw1)
    plt.plot([34.0, 35.0], [0.0, 0.0],    ':', color=lncolor, linewidth=lw1)
    plt.plot([33.5, 35.0], [-0.0575*33.5+0.0901, -0.0575*35.0+0.0901], 
                                          ':', color=lncolor, linewidth=lw1)
    axTS.set_ylim([-2.1,1.5])
    axTS.set_xlim([32.5,34.6])
    axTS.set_ylabel(varlabel[vartitle.index('T')]) 
    axTS.set_xlabel(varlabel[vartitle.index('S')])

    filename = run + '_TS_' + str(z) + m + '_' + placename + '_' + str(startyr) + '-' + str(endyr) 
    if seasonal:
        filename += '_seasonal'
    print(filename)
    plt.savefig(savepath + filename +'.png')
    plt.clf()

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
def tseries1(run_incr,varlist,startyr,endyr,z=0,
             option = 'coord',placename = '',
             latS=latS,lonW=lonW,zab=False,year_overlay=False,
             velocity_vector=False,#runcmp=False,runcmpname='ISMF-noEAIS',
             savepath=savepath_nersc):

    years = np.arange(startyr,endyr,1)
    months = np.arange(1,13,1)
    nt = len(years)*len(months)
    times = np.zeros((nt,))
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run_incr[0])])
    latCell = fmesh.variables['latCell'][:]
    lonCell = fmesh.variables['lonCell'][:]
    xCell = fmesh.variables['xCell'][:]
    yCell = fmesh.variables['yCell'][:]
    if option == 'coord':
        locname = str(latS) + 'S' + str(lonW) + 'W'
        latplt = -1.*latS*deg2rad
        lonplt = (360.-lonW)*deg2rad
        idx = np.argmin( (latCell-latplt)**2 + (lonCell-lonplt)**2)
        placename = ( str(round(latS)) + 'S' + str(round(lonW)) + 'W' ) 
    elif option == 'by_index':
        if placename == 'trough_shelf':
            idx = cells_trough_shelf_lat[-1]-1 
        elif placename == 'trough_ice':
            idx = cells_trough_ice_lat[-1]-1 
    
    # northern limit for subplots
    logical_N = (latCell < lat_N*deg2rad) & (xCell > 0)
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    zice     = fmesh.variables['landIceDraft'][0,:]
    
    if zab:
        m = 'mab'
    else:
        m = 'm'

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
        #if varlim:
        #    plt.clim([varmin[vartitle.index('z')], varmax[vartitle.index('z')]])
        plt.axis('equal')
        cbar = plt.colorbar(cntr1)
        cbar.set_label('Depth (m)')    
        plt.savefig(savepath + 'bathy_' + placename + '.png',dpi=set_dpi)
        plt.close()
    
    data = np.zeros((len(run_incr),len(varlist),nt)) 
    #if runcmp:
    #    data2 = np.zeros((len(varlist),nt)) 
    
    #if 'unormal' in varlist:
    #    _,_,_,transect_angle = pick_transect(option='by_index',
    #                                         run=run,transect_name = 'trough_shelf') 
    t=0
    for yr in years:
        for mo in months:
            times[t] = yr+(mo-1.0)/12.0
       
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            for j,run in enumerate(run_incr):
                filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) 
                           + datestr + '-01.nc')
                f = netCDF4.Dataset(filename, 'r')
                for i,var in enumerate(varlist):
                    if var not in surfvar:
                        zmid,_,_ = zmidfrommesh(fmesh,cellidx=[idx])
                        if zab:
                            zeval = np.add(zmid[0][-1],z)
                        else:
                            zeval = -1*z
                        zidx = np.argmin(np.abs(np.subtract(zmid,zeval)))
 
                    if var in surfvar:
                        data[j,i,t] = f.variables[varname[vartitle.index(var)]][0,idx]
                    else:
                        data[j,i,t] = f.variables[varname[vartitle.index(var)]][0,idx,zidx]
                        #if var == 'unormal':
                f.close()

            #if runcmp:
            #    filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(runcmpname)]) 
            #               + datestr + '-01.nc'
            #    f = netCDF4.Dataset(filename, 'r')
            #    for i,var in enumerate(varlist):
            #        if var in surfvar:
            #            data2[i,t] = f.variables[varname[vartitle.index(var)]][0,idx]
            #        else:
            #            data2[i,t] = f.variables[varname[vartitle.index(var)]][0,idx,zidx]
            #    f.close()    
            t += 1
    
    nrow=len(varlist)
    if velocity_vector:
        nrow += -2
    ncol=1

    fig,axvar = plt.subplots(nrow,ncol,sharex=True)
    t_season.append(1)
    for i in range(len(varlist)):
        if i == nrow-1:
            axvar[i].set(xlabel='year')
        axvar[i].set(ylabel=varlist[i])
        ymin = np.min(data[:,i,:])
        ymax = np.max(data[:,i,:])
        if year_overlay:
            for s,_ in enumerate(season): 
                axvar[i].fill([t_season[s],t_season[s],t_season[s+1],t_season[s+1]],
                              [ymin,ymax,ymax,ymin],
                              facecolor=season_color[s], alpha=0.5, linewidth='none')
    
            for yr in years:
                idx_time = (times>=yr) * (times < yr+1)
                for j,run in run_incr:
                    pc = axvar[i].plot(times[idx_time]-yr,data[j,i,idx_time],
                                       '-', color = run_color[runname.index(run)], alpha = 0.5)
                #if runcmp:
                #    pc2 = axvar[i].plot(times[idx_time]-yr,data2[i,idx_time],'-b',alpha=0.5)
        else:
            for yr in years:
                for s,_ in enumerate(season): 
                    axvar[i].fill([yr + t_season[s],yr + t_season[s],
                                   yr + t_season[s+1], yr + t_season[s+1]],
                                  [ymin,ymax,ymax,ymin],
                                  facecolor=season_color[s], alpha=0.5, linewidth=0)
    
            for j,run in enumerate(run_incr):
                pc = axvar[i].plot(times,data[j,i,:],
                                   '-', 
                                   label = runtitle[runname.index(run)],
                                   color = run_color[runname.index(run)])
            #if runcmp:
            #    pc2 = axvar[i].plot(times,data2[i,:],'-b', label = runname[runname.index(runcmpname)])
    if zab:
        m = 'mab'
    else:
        m = 'm'

    plt.legend(loc=legloc,bbox_to_anchor=bboxanchor)

    #if runcmp:
    #   filename = filename + runcmpname + '_'
    filename = (run_incr[0] + '_' + varlist[0] + varlist[1] + '_t_' + str(z) + m + 
                '_' + placename + '_' + str(startyr) + '-' + str(endyr) ) 
    print(filename)
    plt.savefig(savepath + '/' + filename + '.png',bbox_inches='tight')
    plt.close(fig)

    if velocity_vector:
        i = len(varlist)-2
        plt_aspect = 1/(6*len(years))
        print(plt_aspect)
        width = 10
        fig = plt.figure(figsize=(width,width*plt_aspect*2))
        ax = fig.add_subplot(111)
        ax.set(xlabel='year',ylabel='U (m/s)')
        Umax = np.max(np.sqrt(np.add(np.square(data[i,:]),np.square(data[i+1,:]))))
        y_scalefactor = 1/(12*Umax)
        print(Umax)
        if year_overlay:
            for yr in years:
                idx_time = (times>=yr) * (times < yr+1)
                time = times[idx_time]
                d = data[:,idx_time] 
                if runcmp:
                    d2 = data2[:,idx_time] 
                for ti,t in enumerate(time):
                    plt.plot([t-yr,t-yr+d[i,ti]],[0,d[i+1,ti]],'-k',alpha=0.5)
                    if runcmp:
                        plt.plot([t-yr,t-yr+d2[i,ti]],[0,d2[i+1,ti]],'-b',alpha=0.5)
        else:
            for ti,t in enumerate(times):
                plt.plot([t,t+data[i,ti]],[0,data[i+1,ti]],'-k')
        plt.ylim([-1*Umax,Umax])
        ax.set_aspect(plt_aspect/(2*Umax),adjustable='box')
        #plt.title(run + ': ' + locname)
        filename = run + '_U_t_' + str(z) + m + '_' + placename + '_' + str(startyr) + '-' + str(endyr) 
        print(filename)
        plt.savefig(savepath + '/' + filename + '.png')

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
def hovmoller(run,startyr,endyr,
              option = 'coord', coord=[latS,lonW],
              transect_id = '',
              varlist = ['T','S','rho','u','v'],zlim = [0,-9999],
              limTrue = False,
              input_filename = '',
              savepath = savepath_nersc):

    if option == 'coord':
        years = np.arange(startyr,endyr,1)
        months = np.arange(1,13,1)
        nt = len(years)*len(months)
        times = np.zeros((nt,))
        
        fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
        latCell = fmesh.variables['latCell'][:]
        lonCell = fmesh.variables['lonCell'][:]
        xCell = fmesh.variables['xCell'][:]
        yCell = fmesh.variables['yCell'][:]
        locname = str(round(coord[0])) + 'S' + str(round(coord[1])) + 'W'
        latplt = -1.*coord[0]*deg2rad
        lonplt = (360.-coord[1])*deg2rad
        idx = np.argmin( (latCell-latplt)**2 + (lonCell-lonplt)**2)
        
        # calculate z from depths
        zmid,_,_ = zmidfrommesh(fmesh,cellidx=[idx],vartype='scalar')
        zh = fmesh.variables['layerThickness'][0,idx,:]
        zbottom = np.subtract(zmid[0,:],0.5*zh)
        z = np.zeros((len(zbottom)+1))
        z[0] = zmid[0,0]+zh[0]
        z[1:] = zbottom
        zlim[0] = np.max(z)
        zlim[1] = np.min(z)
        z = z[0:fmesh.variables['maxLevelCell'][idx]]
        nz = len(z)

        data = np.zeros((len(varlist),nz,len(times)))
        
        t=0
        for yr in years:
            for mo in months:
                times[t] = yr+(mo-1.0)/12.0
           
                datestr = '{0:04d}-{1:02d}'.format(yr, mo)
                input_filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) 
                            + datestr + '-01.nc')
                f = netCDF4.Dataset(input_filename, 'r')
                for i,var in enumerate(varlist):
                    data[i,1:,t] = f.variables[varname[vartitle.index(var)]][0,idx,:nz-1]
                
                f.close()
                t += 1
        
        times = np.append(times,np.max(times)+(1/12))
        start_time_idx = 0
        end_time_idx = len(times)+1

    elif option == 'transect':
        #input_filename = run + '_transect_u_'+transect_id+'_' + str(startyr) + '-' + str(endyr)+'.txt'
        if not os.path.exists(savepath + input_filename):
            print(input_filename,' does not exist')
            return
        df = pandas.read_csv(savepath+input_filename)
        #print(df.isnull().sum())
        t = df['decyear'][:]
        times = t.to_numpy(dtype='float32')
        start_time_idx = np.argmin(np.abs(times - startyr))
        end_time_idx = np.argmin(np.abs(times - endyr))
        times = np.append(times,np.max(times)+(1/12))
        # one more time point is needed to specify the right points of quadrilateral
        
        zbottom = np.zeros(200)
        zcol = []
        i = 0
        # ztop is already written to specify the upper points of quadrilateral
        for col in df.columns: 
            if col[0] == '-':
                zcol.append(col)
                zbottom[i] = float(col)
                i += 1
        zbottom = zbottom[zbottom != 0]
        zmid = zbottom[1:] - (zbottom[1:]-zbottom[:-1])/2
        z = zbottom
        nz = len(z)
        data = np.zeros((len(varlist),nz,len(times)))
        for i in range(start_time_idx,end_time_idx):
            data[0,:,i] = df['u_barotropic_sum'][i]
        for i,_ in enumerate(z):
            data[1,i,:-1] = df[zcol[i][:]][:]
        zlim[0] = np.max(z)
        zlim[1] = np.min(z[~np.isnan(data[1,:,0])])
        locname = transect_id
    
    filename = ( run + '_hovmoller_' + varlist[0] + varlist[1] + '_' +
                 locname + '_' + str(startyr) + '-' + str(endyr) )
    nrow=len(varlist)
    ncol=1
    data[np.isnan(data)] = 0
    #print('t = ',np.shape(times[start_time_idx:end_time_idx+1]))
    #print('z = ',np.shape(z))
    #print('data = ',np.shape(data[0,1:,start_time_idx:end_time_idx])) 
    
    fig,axvar = plt.subplots(nrow,ncol,sharex=True)
    for i,var in enumerate(varlist):
        
        cm = plt.get_cmap(varcmap[vartitle.index(var)]) 
        if limTrue:
            cNorm  = colors.Normalize(vmin=varmin[vartitle.index(var)], vmax=varmax[vartitle.index(var)])
        elif var[0] == 'u' or var[0] == 'v':
            vlim = np.max(np.abs(data[i,:,:]))
            cNorm  = colors.Normalize(vmin=-1*vlim, vmax=vlim)
        else:
            cNorm  = colors.Normalize(vmin=np.min(np.abs(data[i,:,:])), 
                                      vmax=np.max(np.abs(data[i,:,:])))
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
 
        pc = axvar[i].pcolormesh(times[start_time_idx:end_time_idx+1], z, 
                                 data[i,1:,start_time_idx:end_time_idx], 
                                 cmap = cm, norm=cNorm) 
        
        if i == nrow-1:
            axvar[i].set(xlabel='Simulated year')
        axvar[i].set(ylabel=varlabel[vartitle.index('z')])
        axvar[i].set_ylim((zlim))
        axvar[i].invert_yaxis()
        
        cbar = fig.colorbar(scalarMap,ax=axvar[i])
        cbar.set_label(varlabel[vartitle.index(var)])
        
        if i == 0:
            print(runtitle[runname.index(run)] + ': ' + loctitle[loc.index(locname)])
            axvar[i].set(title = runtitle[runname.index(run)] + ': ' + loctitle[loc.index(locname)])
    
    print('save plot: ',savepath + filename)
    plt.savefig(savepath + filename + '.png')#,dpi=set_dpi)
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
            maxDepth = -500.,runcmp=False,runcmpname='ISMF-noEAIS',mo = 0,
            savepath=savepath_nersc):

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
    print('idx',idx) 
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
               #f.close()    
               #print(np.shape(data))
               #print(data[0,idx,0])
               #print(z[0])
               if runcmp:
                  filename = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(runcmpname)]) + datestr + '-01.nc'
                  f2 = netCDF4.Dataset(filename, 'r')
                  data2     = f2.variables[varname[vartitle.index(var)]]
                  #f2.close()
               axvar[k].plot(data[0,idx,:], z, 
                             #label="yr{0:04d}".format(yr), 
                             color=c,linestyle=lineStyle[i])
               if runcmp:
                  axvar[k].plot(data2[0,idx,:], z, 
                                #label="yr{0:04d}".format(yr), 
                                color=c,linestyle='--')
    
               axvar[k].set_xlim([varmin[vartitle.index(var)], varmax[vartitle.index(var)]])
               axvar[k].set_ylim([maxDepth,0])
               axvar[k].grid()
    datestr=str(startyr) + '-' + str(endyr)
    plt.title(run + ': ' + datestr)
    filename = savepath + run 
    if runcmp:
       filename = filename + '_' + runcmpname
    filename = filename + '_profiles_' + varlist[0] + varlist[1] + '_' + locname + '_' + datestr
    print(filename)
    plt.savefig(filename + '.png',dpi=set_dpi)
    plt.close()


#----------------------------------------------------------------------
# FLUXGATE
# -- Computes volumetric flux perpendicular to transect line and saves
#    output to a csv file
#
# Inputs:
#   transect_id
#   yrrange 

# optional variables:
#   run       runname, string 
#   runcmp    if true, plot difference between variables of both entries in runname
#   overwrite logical, if true, overwrite image file of same name TODO check
#   plotting  logical, if true, plot map of fluxgate location 
#   savepath  path to save plot images
#----------------------------------------------------------------------
# TODO add option to define fluxgate across topography contour
# TODO add separate function to select transect points, which are then inputs to this function and others
def fluxgate(transect_id, yrrange = [50,51], morange = [1,13],
             run='ISMF', runcmp = False, runcmpname='ISMF-noEAIS',
             overwrite=False, plotting=False, 
             savepath=savepath_nersc):
    
    # import variables from file
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    
    cellidx, idx, dist, transect_angle = pick_transect(option='by_index',
                                    run=run,transect_name = transect_id,overwrite=overwrite) 
                        #lon = [loc_xmin[loc.index(transect_id)],loc_xmax[loc.index(transect_id)]],
                        #lat = [loc_ymin[loc.index(transect_id)],loc_ymax[loc.index(transect_id)]],
    #nLevels = fmesh.dimensions['nVertLevels'][:]
    #cell1 = np.subtract(fmesh.variables['cellsOnEdge'][idx,0],1)
    #cell2 = np.subtract(fmesh.variables['cellsOnEdge'][idx,1],1)
    dv    = fmesh.variables['dvEdge'][idx] # length of edge between vertices
    angle = fmesh.variables['angleEdge'][idx] # angle in rad an edge normal vector makes with eastward
    
    zh   = fmesh.variables['layerThickness'][0,cellidx,:]
    zmax = np.multiply(-1,fmesh.variables['bottomDepth'][cellidx])
    zice = fmesh.variables['landIceDraft'][0,cellidx]
    cell1 = np.subtract(fmesh.variables['cellsOnEdge'][idx,0],1)
    cell2 = np.subtract(fmesh.variables['cellsOnEdge'][idx,1],1)
    #zh1   = fmesh.variables['layerThickness'][0,cell1,:]
    #zh2   = fmesh.variables['layerThickness'][0,cell2,:]
    #zh = (zh1 + zh2)/2
    #zmax1 = fmesh.variables['bottomDepth'][cell1]
    #zmax2   = fmesh.variables['bottomDepth'][cell2]
    #zmax = (zmax1 + zmax2)/-2
    #zice1     = fmesh.variables['landIceDraft'][0,cell1]
    #zice2     = fmesh.variables['landIceDraft'][0,cell2]
    #zice = zice1
    #for i in range(len(cell1)):
    #   zice[i] = np.max([zice1[i],zice2[i]])
    #depths = zmidfrommesh(fmesh,cellidx=idx,vartype='velocity')
    zmid,zbottom,ztop = zmidfrommesh(fmesh,cellidx=cellidx,vartype='scalar')
    #dangle = angle - (transect_angle+(pi/2))
    
    # Create depth vector for reporting cross-transect averaged
    # baroclinic velocities
    zu,idx_unique = np.unique(zbottom,return_index = True)
    temp_idx = np.argsort(-1*zu)
    zu = zu[temp_idx]
    zh_temp = zh.flatten()[idx_unique]
    zuh = zh_temp[temp_idx]
    
    # angle the normal edge velocity makes with transect
    dangle = angle - transect_angle
    dangle[dangle>180] += -360
    edgeSigns = np.divide(dangle,abs(dangle))
 
    #mask = np.zeros(np.shape(depths),dtype=bool)
    #for i,~ in enumerate(mask):
    #    if (depths[i] > zmax[i]) and (depths[i] < zice[i]):
    #        mask[i] = True
    row,col = np.shape(zmid)
    mask = np.zeros((row,col),dtype=bool)
    for i in range(row):
        for j in range(col):
            if (zmid[i,j] > zmax[i]) and (zmid[i,j] < zice[i]):
                mask[i,j] = True
    F = np.multiply(zmid,0)
    
    # compute face area
    area = np.multiply(zmid,0)
    for i in range(row):
        for j in range(col):
            area[i,j] = zh[i,j] * dv[i]

    if plotting: 
        # create mesh variables for plotting
        # distance along transect for plotting
        xpt   = fmesh.variables['xEdge']  [idx]
        ypt   = fmesh.variables['yEdge']  [idx]
        n = np.sqrt( np.square(ypt- ypt[0]) + 
                     np.square(xpt- xpt[0])   )
        yline = np.divide(n,1e3)
        #yfill = np.append(yline[0],yline)
        #yfill = np.append(yfill,yline[-1])
        #sshfill = np.append(0,ssh)
        #sshfill = np.append(sshfill,0)
        #bathymax = np.min(zmax) - 100
        #bathyfill = np.append(bathymax,zmax)
        #bathyfill = np.append(bathyfill,bathymax)
        temp,ymesh= np.meshgrid(np.zeros((col,)),n)
    
    # initialize text files for saving output    
    if overwrite:
        flag='w+'
    else:
        flag='a+'
    
    table_file = open(savepath+run+'_transect_flux_'+transect_id+'_'+str(yrrange[0])+'-'+str(yrrange[1])+'.txt',flag)
    table_file_u = open(savepath+run+'_transect_u_'+transect_id+'_'+str(yrrange[0])+'-'+str(yrrange[1])+'.txt',flag)
    wr = csv.writer(table_file,dialect='excel')
    wru = csv.writer(table_file_u,dialect='excel')
    if runcmp:
        wr.writerow(['year','month','decyear',
                 run+'_flux_pos',run+'_flux_neg',run+'_flux_total',
                 runcmpname+'_flux_pos',runcmpname+'_flux_neg',runcmpname+'_flux_total'])
    else:
        wr.writerow(['year','month','decyear',
                 run+'_flux_pos',run+'_flux_neg',run+'_flux_total'])
   
    col_headings = ['year','month','decyear','F_barotropic','F_baroclinic_pos',
                 'u_barotropic_sum','u_baroclinic_pos',zu[0]+zuh[0]]
    for i in zu:
        col_headings.append(i) 
    wru.writerow(col_headings)
    
    for yr in range(yrrange[0],yrrange[1]):
        for mo in range(morange[0],morange[1]):
            times = (yr+(mo-1.0)/12.0)

            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) 
                       + datestr + '-01.nc')
            print(run,': ',filename)
            
            if not os.path.exists(filename):
              print('does not exist')
              continue

            f = netCDF4.Dataset(filename, 'r')
            ssh   = f.variables['timeMonthly_avg_ssh'][0,cellidx]
            if 'timeMonthly_avg_normalTransportVelocity' in f.variables.keys():
              u= (f.variables['timeMonthly_avg_normalTransportVelocity'][0,idx,:])
            elif 'timeMonthly_avg_normalVelocity' in f.variables.keys():
              u= (f.variables['timeMonthly_avg_normalVelocity'][0,idx,:])
              if 'timeMonthly_avg_normalGMBolusVelocity' in f.variables.keys():
                u += (f.variables['timeMonthly_avg_normalGMBolusVelocity'][0,idx,:])
            else:
              raise KeyError('no appropriate normalVelocity variable found')
            for i,_ in enumerate(idx):
                u[i,:] = edgeSigns[i]*u[i,:]
            F = np.multiply(area,u)
            
            if runcmp:
                filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(runcmpname)])
                            + datestr + '-01.nc')
                f2 = netCDF4.Dataset(filename, 'r')
                u2 = f2.variables['timeMonthly_avg_normalVelocity'][0,idx,:]
                for i,_ in enumerate(idx):
                    u2[i,:] = edgeSigns[i]*u2[i,:]
                F2 = np.multiply(np.multiply(area,u2),m3ps_to_Sv)
            else:
                F2 = [0]
            
            # Compute barotropic/baroclinic velocities column-wise
            u_barotropic = np.zeros(np.shape(u))
            u_baroclinic= np.zeros(np.shape(u))
            column_flux = np.zeros(np.shape(cellidx))
            column_area = np.zeros(np.shape(cellidx))
            for i,_ in enumerate(idx):
                column_flux[i] = np.sum(F[i,:]*mask[i,:])
                column_area[i] = np.sum(np.multiply(area[i,:],mask[i,:]))
                u_barotropic[i,:] = np.divide(column_flux[i],
                                              np.add(column_area[i],1e-10))
                u_baroclinic[i,:] = np.subtract(u[i,:],
                                                u_barotropic[i,:])
            u_barotropic = np.multiply(u_barotropic,mask)
            u_baroclinic = np.multiply(u_baroclinic,mask)
            F_barotropic = np.multiply(u_barotropic,area) 
            F_baroclinic = np.multiply(u_baroclinic,area)
            
            # check that F_baroclinic_sum = 0
            for i,_ in enumerate(idx):
                if abs(np.sum(F_baroclinic[i,:])) > 1e-10:
                    print('Column-integrated baroclinic flux exceeds threshold')
            
            # compute cross-transect averaged barotropic and baroclinic velocities
            u_baroclinic_zsum = np.zeros((np.shape(zu)))
            width_zsum = np.zeros((np.shape(zu)))
            for i,_ in enumerate(cellidx):
                for j,_ in enumerate(zmid):
                    for k,_ in enumerate(zu):
                        if (zbottom[i,j] >= zu[k]) and (zbottom[i,j] < zu[k]+zuh[k]):
                            u_baroclinic_zsum[k] += u_baroclinic[i,j] * dv[i]
                            width_zsum[k] += dv[i]
            u_baroclinic_zsum = np.divide(u_baroclinic_zsum,width_zsum)

            F_barotropic_sum = np.sum(F_barotropic)
            F_baroclinic_sum = np.sum(F_baroclinic)
            F_baroclinic_mag = np.sqrt(np.sum(np.square(F_baroclinic)))
            F_baroclinic_pos = np.sum(F_baroclinic[F_baroclinic>0])
            area_sum = np.sum(area)
            # the depth-averaged velocities
            u_barotropic_sum = F_barotropic_sum/area_sum
            u_baroclinic_sum = F_baroclinic_pos/area_sum # TODO change to mag or pos
            print('F_barotropic_sum = ',F_barotropic_sum)
            print('F_baroclinic_sum = ',F_baroclinic_sum)
            print('F_baroclinic_mag = ',F_baroclinic_mag)
            print('F_baroclinic_pos = ',F_baroclinic_pos)
            print('u_barotropic_sum = ',u_barotropic_sum)
            print('u_baroclinic_zsum = ',u_baroclinic_zsum)
            row_entries = [yr,mo,times,
                           F_barotropic_sum * m3ps_to_Sv,
                           F_baroclinic_pos * m3ps_to_Sv,
                           u_barotropic_sum, u_baroclinic_sum,
                           0]
            for i in u_baroclinic_zsum:
                row_entries.append(i)
            wru.writerow(row_entries)
            
            F = F[mask]
            Fpos = np.sum(F[F>0])
            Fneg = np.sum(F[F<0])
            Fsum = np.sum(F)
            if runcmp:
                F2 = F2[mask]
                F2pos = np.sum(F2[F2>0])
                F2neg = np.sum(F2[F2<0])
                F2sum = np.sum(F2)
            else:
                F2pos = 0
                F2neg = 0
                F2sum = 0
            if runcmp:
                wr.writerow([yr,mo,times,
                             Fpos*m3ps_to_Sv,
                             Fneg*m3ps_to_Sv,
                             Fsum*m3ps_to_Sv,
                             F2pos*m3ps_to_Sv,
                             F2neg*m3ps_to_Sv,
                             F2sum*m3ps_to_Sv])
            else:
                wr.writerow([yr,mo,times,
                             Fpos*m3ps_to_Sv,
                             Fneg*m3ps_to_Sv,
                             Fsum*m3ps_to_Sv])
            
            # in a separate text file, write 1st row depths, subsequent rows per time
            
            if plotting: 
                ssh      = f.variables['timeMonthly_avg_ssh'][0,cell1]
                fig2 = plt.figure()
                filename = run
                if runcmp:
                    filename = filename + '_cmp'
                filename = filename + '_u_' + transect_id + '_' + datestr
                plt.scatter(np.divide(ymesh[mask],1e3),
                            -1*(depths[mask]), 
                            s = 1,c=u[mask],
                            cmap = "cmo.balance",vmin=-1*np.max(abs(u[mask])),vmax=np.max(abs(u[mask])))
                plt.plot(yline, -1*zice, 
                         color = 'blue', marker = '.', linestyle = '-')
                plt.plot(yline, -1*zmax, 
                         color = 'black', marker = '.', linestyle = '-')
                ax = plt.gca()
                ax.invert_yaxis()
                cbar = plt.colorbar()
                cbar.set_label('U')
                cbar.set_clim(-1*np.max(abs(u[mask])),np.max(abs(u[mask])))
                plt.xlabel('Distance (km)')
                plt.ylabel('Depth (m)')
                plt.title(run + ': ' + datestr)
                plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
                print(filename) 
                plt.close()
                fig2 = plt.figure()
                filename = run
                if runcmp:
                    filename = filename + '_cmp'
                filename = filename + '_f_' + transect_id + '_' + datestr
                plt.scatter(np.divide(ymesh[mask],1e3),
                            -1*(depths[mask]), 
                            s = 1,c=F[mask],
                            cmap = "cmo.balance",vmin=-1*np.max(abs(F[mask])),vmax=np.max(abs(F[mask])))
                plt.plot(yline, -1*zice, 
                         color = 'blue', marker = '.', linestyle = 'none')
                plt.plot(yline, -1*zmax, 
                         color = 'black', marker = '.', linestyle = 'none')
                ax = plt.gca()
                ax.invert_yaxis()
                cbar = plt.colorbar()
                cbar.set_label('F/a')
                #cbar.set_clim(-1*np.max(abs(F[mask])),np.max(abs(F[mask])))
                plt.xlabel('Distance (km)')
                plt.ylabel('Depth (m)')
                plt.title(run + ': ' + datestr)
                plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
                print(filename) 
                plt.close()


#----------------------------------------------------------------------
# TRANSECT
# -- plot variable in color vs. depth across a profile in lon or lat
# -- if runcmp, plot difference in two fields
#
# Inputs:
#   pick_option  'by_index'
#   latS      latitude range, always in Southern Hem, positive, 
#             vector of length 2, real 
#   lonW      longitude range, always in Western Hem, positive,
#             vector of length 2, real 
#   varlist   variables to plot, list of strings
#   yr        lower limit on simulated year to plot, real
#   mo        month to plot data from, if 0 then plot all months of data

# optional variables:
#   varlim
#   zscale    scale of depth axis
#   run       runname, string 
#   runcmp    if true, plot difference between variables of both entries in runname
#   new       
#   ops       operations to perform on variables in varlist
#   savepath  path to save plot images
#----------------------------------------------------------------------
def transect(pick_option, yr_incr, mo_incr, varlist, 
             transect_name = '',plot_transect=False,
             plot_method = 'tricontourf',
             lat=[latmin,latmax],lon=[lonmin,lonmax], 
             varlim = False, normal = False, zscale = 'linear', 
             run='ISMF', runcmp = False, runcmpname = 'ISMF-noEAIS',
             overwrite = False, ops = [''],
             zlim = [-9999,-9999],
             save_transect_mean = False,
             savepath=savepath_nersc,figure_format = 'png'):
    
    if ops[0] == '':
       ops = ['' for i in varlist]
    
    if pick_option == 'coord':
        transect_name = (str(int(abs(lat[0]))) + 'S' + 
                         str(int(abs(lon[0]-360))) + 'W-' + 
                         str(int(abs(lat[1]))) + 'S' + 
                         str(int(abs(lon[1]-360))) + 'W' )
    
    # create mesh variables for plotting
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    cellidx,edgeidx,dist,angle = pick_transect(option = pick_option,
                                       lat = lat, lon = lon, 
                                       transect_name = transect_name,
                                       overwrite = plot_transect) 
    dist = np.divide(dist,1e3) 
    # constants
    #dlat = 0.15 # at 30km resolution, distance between cells in latitude space
    #dlon = 0.98
    
    # import variables from file
    #latCell  = fmesh.variables['latCell'][idx]
    #lonCell  = fmesh.variables['lonCell'][idx]
    #idxCell  = fmesh.variables['indexToCellID'][idx]
    #xCell    = fmesh.variables['xCell'][:]
    #yCell    = fmesh.variables['yCell'][:]
    kmax     = fmesh.variables['maxLevelCell'][cellidx]
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][cellidx])
    zh       = fmesh.variables['layerThickness'][0,cellidx]
    #icemask  = fmesh.variables['landIceMask'][idx]
    zice     = fmesh.variables['landIceDraft'][0,cellidx]
    bathymax = np.min(zmax) - 100
    
    # calculate z from depths
    zmid,zbottom,ztop = zmidfrommesh(fmesh, cellidx = cellidx, 
                         vartype = 'scalar')
    zmesh = np.zeros((len(cellidx),len(zbottom[0,:])+2))
    if plot_method == 'tricontourf':
        zmesh[:,1:-1] = zmid[:,:]
        for idx,_ in enumerate(cellidx):
           zmesh[idx,kmax[idx]+1:] = zbottom[idx,kmax[idx]]
        zmesh[:,0] = ztop[:,0]
    if plot_method == 'pcolormesh':
        zmesh[:,0] = ztop[:,0]
        zmesh[:,1:] = zbottom[:,1:]

    _,ymesh  = np.meshgrid(zmesh[0,:], dist)
    #print('zmesh shape = ',np.shape(ymesh))
    #print('ymesh shape = ',np.shape(ymesh))
    #print('ymesh[0] = ',ymesh[0,:])

    # TODO only used for plotting, ignore if not plotting 
    yfill = np.append(dist[0],dist)
    yfill = np.append(yfill,dist[-1])
    bathyfill = np.append(bathymax,zmax)
    bathyfill = np.append(bathyfill,bathymax)

    for yr in yr_incr:
        for mo in mo_incr:
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'
                        .format(runpath[runname.index(run)]) 
                        + datestr + '-01.nc')
            f = netCDF4.Dataset(filename, 'r')
            if runcmp:
                filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'
                            .format(runpath[runname.index(runcmpname)]) 
                            + datestr + '-01.nc')
                f2 = netCDF4.Dataset(filename, 'r')

            
            ssh      = f.variables['timeMonthly_avg_ssh'][0,cellidx]
            sshfill = np.append(0,np.minimum(ssh,zice))
            sshfill = np.append(sshfill,0)
            
            # TODO finish this section
            #ddist = np.zeros((len(dist)+1))
            #ddist[0] = ddist[0] - (dist[1]-dist[0])
            #ddist[1:-1] = dist[1:] - (dist[1:]-dist[:-1])
            #ddist[-1] = dist[-1] + (dist[-1]-dist[-2])
            #print('ddist shape = ',np.shape(ddist))
            #zzmesh = np.zeros((len(cellidx),len(zh[0,:]+1)))
            #print('zmesh shape = ',np.shape(zmesh))
            #print('zh shape = ',np.shape(zh))
            #print('zzmesh shape = ',np.shape(zzmesh))
            #zzmesh[:,:-1] = zmesh - np.multiply(0.5,zh)
            #zzmesh[:,-1]  = ssh
            #print(zzmesh[0,:])
    
            for var in varlist:
                image_filename = savepath + run
                if runcmp:
                    image_filename = image_filename + '_cmp'
                image_filename = image_filename + '_' + var 
                if var == 'u' and normal :
                   image_filename = image_filename + 'normal' 
                image_filename = ( image_filename + 
                                   ops[varlist.index(var)] +  
                                   '_' + transect_name + '_' + datestr + '_lim' + 
                                   str(varlim)) 

                if not overwrite and os.path.exists(image_filename + '.' + figure_format):
                    print(image_filename + ' exists')
                    continue
                
                data_import = f.variables[varname[vartitle.index(var)]][0,cellidx,:]
                # convert velocity to the normal velocity
                if var == 'u' and normal:
                    u = data_import
                    v = f.variables[varname[vartitle.index('v')]][0,cellidx,:]
                    u_angle = np.arctan2(v,u)
                    transect_angle_normal = np.add(angle, pi/2)
                    u_norm = np.sqrt(np.add(np.square(u),np.square(v)))
                    data_import = np.multiply(u_norm, 
                                       np.cos(np.subtract(u_angle,
                                              transect_angle_normal)))
                
                #TODO, test
                if runcmp: 
                    if var == 'u' and normal:
                        print('normal velocity plot not yet compatible with difference between runs')
                        return
                    data2 = f2.variables[varname[vartitle.index(var)]][0,cellidx,:]
                    data_import = np.subtract(data_import,data2)
                
                data = np.zeros((np.shape(zmesh)))
                if plot_method == 'tricontourf':
                    data[:,1:-1] = data_import
                    data[:,-1] = data[:,-2]
                    data[:,0] = data[:,1]
                    # for interpolation, add the variable assigned at bottomDepth
                    for idx,i in enumerate(cellidx):
                       data[idx,kmax[idx]+1] = data_import[idx,kmax[idx]]
                if plot_method == 'pcolormesh':
                    data[:,1:] = data_import
                    data[:,0] = data[:,1]
                
                mask = np.zeros((np.shape(zmesh)),dtype=bool)
                mask = (data < bad_data)
                for idx,i in enumerate(cellidx):
                    mask[idx,:] = (zmesh[idx,:] > zice[idx]) | (zmesh[idx,:] < zmax[idx])
                    mask[idx,kmax[idx]+1:] = True
                data_masked  = data[~mask]
                zmesh_masked = zmesh[~mask]
                ymesh_masked = ymesh[~mask]
               
                # TODO revise and test 
                if ops[varlist.index(var)] == 'barotropic':
                    for idx,i in enumerate(cellidx):
                        u_depth_mean = 0
                        for jdx,j in enumerate(depths):
                            #wct = 
                            u_depth_mean += (layer_thickness/wct) * data_masked[idx,jdx]
                            data_masked[idx,:] = np.mean(temp[~mask[idx,:]])
                    
                elif ops[varlist.index(var)] == 'baroclinic':
                    for idx,i in enumerate(cellidx):
                        temp = data_masked[idx,:]
                        data_masked[idx,:] = np.subtract(temp,
                                                np.mean(temp[~mask[idx,:]]))
                
                # plots
                
                clevels = np.arange(varmin[vartitle.index(var)], 
                                    varmax[vartitle.index(var)],
                                    dvar[vartitle.index(var)]   )
                if clevels[0] > np.min(data_masked.flatten()):
                    clevels = np.append(np.min(data_masked.flatten()),clevels)
                if clevels[-1] < np.max(data_masked.flatten()):
                    clevels = np.append(clevels,np.max(data_masked.flatten()))
                
                fig = plt.figure()
                
                if runcmp:
                    cmap1 = "cmo.balance"
                else:
                    cmap1 = varcmap[vartitle.index(var)]
                
                if plot_method == 'tricontourf':
                    cntr2 = plt.tricontourf(np.transpose(ymesh_masked.flatten()), 
                                        np.transpose(
                                           np.abs(zmesh_masked.flatten())), 
                                        np.transpose(
                                           data_masked.flatten()), 
                                        levels=clevels,cmap=cmap1)
                # TODO pcolormesh not ready yet
                elif plot_method == 'pcolormesh':
                    pc = plt.pcolormesh(yymesh, zzmesh, 
                                        data_mesh,
                                        cmap = cmap1, norm=cNorm) 
                plt.plot(dist, np.abs(ssh), 
                         color = 'black', marker = '.', linestyle = '-')
                plt.plot(ymesh[~mask].flatten(), 
                         np.abs(zmesh[~mask].flatten()), 
                         '.', color = 'white', markersize = 1)#, fillstyle = 'none')
                plt.fill(yfill, np.abs(sshfill), c = 'white', alpha = 1)
                plt.plot(dist, np.abs(zmax), 
                         color = 'black', marker = '.', linestyle = '-')
                plt.fill(yfill, np.abs(bathyfill), c = 'grey', alpha = 1)
                ax = plt.gca()

                if varlim:
                    if runcmp:
                        if var == 'T':
                            cmax = 1.
                        else:
                            cmax = max(abs(np.percentile(data_masked[~mask].flatten(),
                                                         10)),
                                       abs(np.percentile(data_masked[~mask].flatten(),
                                                         90)))
                        plt.clim([-1*cmax,cmax])
                    else:
                        plt.clim([varmin[vartitle.index(var)], 
                                  varmax[vartitle.index(var)]])

                if zscale == 'log':
                    ax.set_yscale('log')
                if zlim[0] != -9999:
                    plt.ylim(np.multiply(zlim,-1))
                    image_filename += '_zmax' + str(int(-1*zlim[1]))
                else:
                    plt.ylim([np.max(sshfill),np.abs(np.min(zmax))])
                ax.invert_yaxis()
                cbar = plt.colorbar()
                cbar.set_label(var)
                plt.xlabel('Distance (km)')
                plt.ylabel('Depth (m)')
                #plt.xlim([np.min(dist)/1e3,np.max(dist)/1e3])
                plt.title(runtitle[runname.index(run)] + ': ' + datestr)
                
                plt.savefig(image_filename+'.'+figure_format,dpi=set_dpi)
                print(image_filename) 
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
def plot_zice_map(run = 'ISMF',locname = 'fris',savepath=savepath_nersc):

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

def calc_stresscurl_t(var,run=['ISMF'],yrrange=[70,71],
                      locname='wedwang',coord='lat',
                      overwrite=False,map_output=False,
                      savepath=savepath_nersc):
    filename = run[0]
    if len(run) > 1:
        filename = filename + '_cmp_' + run[1]
    filename = filename + '_windstresscurl_' + locname 
    print(filename)
    
    if locname not in loc:
        print('locname is not defined')
        return
    
    # check if plot was already generated
    if os.path.exists(savepath + filename + '.png'):
        print('file exists')
        if not overwrite:
            print('skipping file')
            return
    
    xmin = loc_xmin[loc.index(locname)]
    xmax = loc_xmax[loc.index(locname)]
    ymin = loc_ymin[loc.index(locname)]
    ymax = loc_ymax[loc.index(locname)]
    # open data files
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    # import variables from file
    xCell    = fmesh.variables['xCell'][:]
    yCell    = fmesh.variables['yCell'][:]
    latCell  = fmesh.variables['latCell'][:]
    lonCell  = fmesh.variables['lonCell'][:]
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    zice     = fmesh.variables['landIceDraft'][0,:]
    
    idx_A = np.argwhere((xCell < 2.5e6) & (xCell > 0) & (yCell < 0) & (yCell > -1.5e6) & (latCell < lat_N*deg2rad))
    logical_A = idx_A[:,0]
    print('ndom=',len(logical_A)) 
    if coord == 'xy':
        idx_loc = np.argwhere((xCell < xmax) & (xCell > xmin) 
                               & (yCell < ymax) & (yCell > ymin))
    elif coord == 'lat':
        idx_loc = np.argwhere((lonCell < xmax*deg2rad) & (lonCell > xmin*deg2rad) 
                           & (latCell < ymax*deg2rad) & (latCell > ymin*deg2rad))
    
    logical_N = idx_loc[:,0]
    x = xCell[logical_N]
    y = yCell[logical_N]
    print('nx = ',len(x),', ny = ',len(y))
    
    dx = (np.max(x)-np.min(x))/sqrt(len(x))
    dy = (np.max(y)-np.min(y))/sqrt(len(x))
    #dx = np.min(x - x[0])
    #dy = np.min(y - y[0])
    print('dx= ',dx,', dy = ',dy)
    x_reg = np.arange(np.min(x),np.max(x)+dx,dx)
    y_reg = np.arange(np.min(y),np.max(y)+dx,dy)
    print('x2 = ',len(x_reg),', y2 = ',len(y_reg))
    xi,yi = np.meshgrid(x_reg,y_reg)
    print('xi = ',np.shape(xi))
    #reg_vector_x = np.ndarray.flatten(reg_mesh_x)
    #reg_vector_y = np.ndarray.flatten(reg_mesh_y)
    #print('nx2 = ',len(reg_vector_x),', ny2 = ',len(reg_vector_y))

    if not os.path.exists(savepath + 'bathy_' + locname + '.png'):
        fig = plt.figure()
        cntr1 = plt.tricontourf(yCell[logical_A].flatten(), xCell[logical_A].flatten(), 
                                zmax[logical_A].flatten(), 20, cmap="viridis")
        plt.plot(y,x, 'o', color = 'white', 
                 markersize = 4, fillstyle = 'none')#, alpha = 0.5)
        plt.plot(yi,xi, 'o', color = 'red', 
                 markersize = 4)
        cntr = plt.tricontour(yCell[logical_A].flatten(), xCell[logical_A].flatten(), 
                               zice[logical_A].flatten(), [-10], colors = 'k')
        plt.clim([varmin[vartitle.index('z')], varmax[vartitle.index('z')]])
        plt.axis('equal')
        cbar = plt.colorbar(cntr1)
        cbar.set_label('Depth (m)')    
        plt.savefig(savepath + 'bathy_' + locname + '.png',dpi=set_dpi)
        plt.close()
    curl_t=[]
    if runcmp:
        curl_t_cmp=[]
    if overwrite:
        flag='w+'
    else:
        flag='a+'
    
    table_file = open(savepath+'windstresscurl_'+locname+'_'+str(yrrange[0])+'-'+str(yrrange[1])+'.txt',flag)
    wr = csv.writer(table_file,dialect='excel')
    headings = ['year','month','decyear']
    for run_incr in run:
        headings.append(i+'_curl')
    wr.writerow(headings)
    curl = np.zeros((3,))
    for yr in range(yrrange[0],yrrange[1]):
        for mo in range(1,13):
            times = (yr+(mo-1.0)/12.0)
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            for k,run_incr in enumerate(run):
                filein = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run_incr)]) 
                          + datestr + '-01.nc')
                f = netCDF4.Dataset(filein, 'r')
                taux = f.variables['timeMonthly_avg_windStressZonal'][0,logical_N]
                tauy = f.variables['timeMonthly_avg_windStressMeridional'][0,logical_N]
                taux_i = interp.griddata((x,y),taux,(xi,yi),method='linear')
                tauy_i = interp.griddata((x,y),tauy,(xi,yi),method='linear')
                dtauy_dx = nan*np.ones((len(x_reg)-1,len(y_reg)-1))
                dtaux_dy = nan*np.ones((len(x_reg)-1,len(y_reg)-1))
                for i in range(1,len(x_reg)):
                    dtauy_dx[i-1,:] = np.divide(tauy_i[i,1:] - tauy_i[i-1,1:],dy)
                for j in range(1,len(y_reg)):
                    dtaux_dy[:,j-1] = np.divide(taux_i[1:,j] - taux_i[1:,j-1],dx)
                curl[k] = np.nanmean(np.subtract(dtauy_dx,dtaux_dy))
            wr.writerow([yr,mo,times,curl])

            if map_output:
                ms = 50
                fig = plt.figure()
                ax1 = fig.add_subplot(121)
                ax1.scatter(x,y,s=ms,c=tauy,edgecolors='k')
                ax1.scatter(xi,yi,s=ms,c=tauy_i,edgecolors='r')
                ax1.set_xlabel('xi',fontsize=16)
                ax1.set_ylabel('yi',fontsize=16)
                
                ax2 = fig.add_subplot(122)
                ax2.scatter(x,y,s=ms,c=taux,edgecolors='k')
                ax2.scatter(xi,yi,s=ms,c=taux_i,edgecolors='r')
                ax2.set_xlabel('xi',fontsize=16)
                ax2.set_ylabel('yi',fontsize=16)
                plt.savefig(savepath+run_incr+'tau_xy_sc.png',dpi=100)
                plt.close(fig)
        
    return
# filename list of filenames corresponding to runs
# run      list of runnames of same length as filename
def plot_stresscurl_t(filename,run = ['ISMF'], tlim=[9999.,9999.],
                      plot_difference = False, savepath=savepath_nersc):
    plot_filename = ''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for k, run_incr in enumerate(run):
        df = pandas.read_csv(savepath + filename[k] + '.txt')
        t = df['decyear'][:]
        curl = df[run_incr+'_curl'][:]
        if not plot_difference:
            ax.plot(t,curl,
                    color=run_color[runname.index(run_incr)],
                    label=runtitle[runname.index(run_incr)])
        if plot_difference and k == 0:
            curl_cmp = curl
        plot_filename = plot_filename + run_incr + '_'

    if plot_difference:
        dcurl = np.subtract(curl_cmp,curl)
        ax.plot(t,dcurl,label = runname[0]+'-'+runname[1])
        plot_filename = plot_filename + 'cmp'

    ax.set_xlabel(r'Year',fontsize=fs)
    if tlim[0] != 9999.:
       ax.set_xlim(tlim)
       plot_filename = plot_filename + str(int(tlim[0])) +'-'+ str(int(tlim[1]))

    ax.set_ylabel(r'Wind stress curl (N m$^{-3}$)',fontsize=fs)

    ax.legend(loc=9,bbox_to_anchor=(0.15,-0.15),fontsize=fs)
    plt.savefig(savepath+filename+'.png',dpi=100,bbox_inches='tight')
    plt.close(fig)
    return

#------------------------------------------------------------------------------
# PLOT_FLUXGATE_T 
# -- Opens a csv file with mean volumetric flux through the fluxgate as 
#    a function of time and generates a line plot
# 
# Inputs:
# 
#------------------------------------------------------------------------------
def plot_fluxgate_t(filename,tlim=[9999.,9999.],run_incr=['ISMF'],runcmpname='ISMF-noEAIS',
                    mode = 'barotropic-baroclinic',savepath=savepath_nersc):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if mode == 'pos-neg':
        #print(filename)
        df = pandas.read_csv(savepath+filename+'.txt')
        #for col in df.columns: 
        #   print(col) 
        print(df.isnull().sum())
        t = df['decyear'][:]
        for run in run_incr:
            Fpos = df[run+'_flux_pos'][:]          
            Fneg = df[run+'_flux_neg'][:]          
            Fsum = df[run+'_flux_total'][:]        
            F2pos = df[runcmpname+'_flux_pos'][:]  
            F2neg = df[runcmpname+'_flux_neg'][:]  
            F2sum = df[runcmpname+'_flux_total'][:]
            ax.plot(t, Fpos,'--b',label=run+' +')
            ax.plot(t, Fneg, ':b',label=run+' -')
            ax.plot(t, Fsum, '-b',label=run+' total')
            ax.plot(t,F2pos,'--k',label=runcmpname+' +')
            ax.plot(t,F2neg, ':k',label=runcmpname+' -')
            ax.plot(t,F2sum, '-k',label=runcmpname+' total')
    elif mode == 'barotropic-baroclinic':
        t_season.append(1)
        barotropic_ls = '-'
        baroclinic_ls = '--'
        ylim = [0,0]
        for run in run_incr:
            input_filename = run + filename[4:] + '.txt'
            df = pandas.read_csv(savepath+input_filename)
            t = df['decyear'][:]
            F_barotropic = df['F_barotropic'][:]          
            F_baroclinic = df['F_baroclinic_pos'][:]          
            ax.plot(t, F_barotropic, 
                    linestyle = barotropic_ls, color = run_color[run_incr.index(run)],
                    label = runtitle[runname.index(run)] + ',barotropic')
            ax.plot(t, F_baroclinic, 
                    linestyle = baroclinic_ls, color = run_color[run_incr.index(run)],
                    label = runtitle[runname.index(run)] + ',baroclinic')
            ylim = [np.min([ylim[0],np.min(F_barotropic),np.min(F_baroclinic)]),
                    np.max([ylim[1],np.max(F_barotropic),np.max(F_baroclinic)])]
        for yr in range(floor(np.min(t)),ceil(np.max(t))):
            for s,_ in enumerate(season): 
                ax.fill([yr + t_season[s],   yr + t_season[s],
                         yr + t_season[s+1], yr + t_season[s+1]],
                         [ylim[0],ylim[1],ylim[1],ylim[0]],
                         facecolor=season_color[s], alpha=0.5, linewidth=0)
    
    if tlim[0] != 9999.:
       ax.set_xlim(tlim)
       filename = filename + '_' + str(int(tlim[0])) +'-'+ str(int(tlim[1])) + '_season'
    ax.set_xlabel(r'Year',fontsize=fs)
    ax.set_ylabel(r'Volume Flux (Sv)',fontsize=fs)
    #ax.set_ylabel(r'Volume Flux (m$^{3}$s$^{-1})',fontsize=fs)
    ax.legend(loc=9,bbox_to_anchor=(0.15,-0.15),fontsize=fs)
    plt.savefig(savepath+filename+'.png',dpi=100,bbox_inches='tight')
    plt.close(fig)
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
def plot_surf_var(var,yr,mo,run=['ISMF'],locname='fris',plottype = 'abs',
                  z=0,zab=False,level=bad_data,
                  varlim=False,overwrite=False,
                  savepath=savepath_nersc):

    if len(run) == 2:
        runcmp = True
    else:
        runcmp = False
    if locname not in loc:
        print('locname is not defined')
        return
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
    filename = run[0]
    if runcmp:
        filename = filename + '_cmp_' + run[1]
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
    fmesh = netCDF4.Dataset(meshpath[runname.index(run[0])])
    filein = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run[0])]) + datestr + '-01.nc'
    f = netCDF4.Dataset(filein, 'r')
    if runcmp:
       filein2 = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run[1])]) + datestr + '-01.nc'
       f2 = netCDF4.Dataset(filein2, 'r')
    
    # import variables from file
    latCell  = fmesh.variables['latCell'][:]
    lonCell  = fmesh.variables['lonCell'][:]
    idxCell  = fmesh.variables['indexToCellID'][:]
    xCell    = np.divide(fmesh.variables['xCell'][:],1e3)
    yCell    = np.divide(fmesh.variables['yCell'][:],1e3)
    kmax     = fmesh.variables['maxLevelCell'][:]
    zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    zh       = fmesh.variables['layerThickness'][0,:]
    zice     = fmesh.variables['landIceDraft'][0,:]
    landicemask  = fmesh.variables['landIceMask'][0,:]
    
    # northern limit for subplots
    idx_N = np.argwhere((xCell<loc_xmax[loc.index(locname)]/1e3) & 
                        (xCell > 0) & 
                        (yCell < 0) & 
                        (yCell > loc_ymin[loc.index(locname)]/1e3) & 
                        (latCell < lat_N*deg2rad))
    logical_N = idx_N[:,0]
    x = xCell[logical_N]
    y = yCell[logical_N]
    idx = idxCell[logical_N]
    
    # calculate z from depths
    zmid,_,_ = zmidfrommesh(fmesh,cellidx=idx)
    
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
                    #zidx[i] = int(np.argmin(np.abs(np.subtract(zmid[idx,:],zeval))))
                    zidx[i] = int(np.argmin(np.abs(np.subtract(zmid[i,:],zeval))))
        else:
            for i,idx in enumerate(logical_N):
                #zidx[i] = int(np.argmin(np.abs(np.subtract(zmid[idx,:],-1*z))))
                zidx[i] = int(np.argmin(np.abs(np.subtract(zmid[i,:],-1*z))))
        
        # get data
        icemask = (landicemask[logical_N] == 1)
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
            bad_idx = ((dataz > bad_data) | (landicemask[logical_N] == 1))
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
            if var == 'ssh' and runcmp:
                plt.clim([varmin[vartitle.index(var+'_cmp')], varmax[vartitle.index(var+'_cmp')]])
        elif level > bad_data and not runcmp:
            #plt.clim([0,2000])
            plt.clim([np.percentile(dataz[bad_idx],10),np.percentile(dataz[bad_idx],90)])
        #else:
        #    plt.clim([np.percentile(dataz[bad_idx],10),np.percentile(dataz[bad_idx],90)])
        
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
        #gl1 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
        #                     zice[logical_N].flatten(), [-1], colors = 'b', linewidths = 2)
        #gl1 = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
        #                     landicemask[logical_N].flatten(), [0], colors = 'b', linewidths = 2)
        
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
            cbar1.set_label(varlabel[vartitle.index(var)])
        fig.tight_layout()
    
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
    
    if runcmp:
        plt.title(run[0] + ' - ' +run[1] + ': ' + datestr)
    else:
        plt.title(run[0] + ': ' + datestr)
    
    #ax.set_ylim([loc_ymin[loc.index(locname)]/1e3,loc_ymax[loc.index(locname)]/1e3])
    #ax.set_xlim([loc_xmin[loc.index(locname)]/1e3,loc_xmax[loc.index(locname)]/1e3])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.savefig(savepath + '/' + filename + '.png',dpi=set_dpi)
    plt.close()

