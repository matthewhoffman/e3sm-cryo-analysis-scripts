#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

import sys
import os
import csv
import netCDF4
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
from math import pi,sqrt,nan,atan,cos,floor,ceil
import pandas

global bad_data, bad_data2, deg2rad, lat_N, runname, runpath, meshpath, vartitle, varname, varmin, varmax, varcmap, surfvar, dvar

pltlib.rc_file('rcparams.txt', use_default_template=True)

m3ps_to_Sv = 1e-6 # m^3/sec flux to Sverdrups
bad_data = -1e33
bad_data2 = 0.
deg2rad  = pi/180.
lat_N = -50 # northern limit of domain
    
loc = ['frisEAcoast','fris','wedwang','trough_shelf','trough_ice','76S31W','75S31W']
cells_trough_ice_lat = [33000,77033,11788,129098,49348] # index is value - 1, located just under ice
cells_trough_shelf_lat = [165785,166563,130260,191987,85569] # index is value - 1, located just under ice
#cells_trough_shelf_lat = [224105,165785,166563,130260,191987] # index is value - 1, located just under ice
cells_trough_shelf2_lat = [152383,144943,67352] # index is value - 1, located just under ice
trans_opt = ['' for i in loc]
#trans_opt[loc.index('trough_shelf')] = 'coord'
#trans_opt[loc.index('trough_ice')] = 'coord'
loctitle = ['' for i in loc]
loctitle[loc.index('trough_shelf')] = 'trough, shelf break'
loctitle[loc.index('trough_ice')] = 'trough, under ice shelf front'
loctitle[loc.index('76S31W')] = 'M31W'
loctitle[loc.index('75S31W')] = 'S4E'

# location M31W from Ryan et al. 2017 JGR Oceans
latS = 76 + 2.75/60
lonW = 30 + 59.65/60

latmin = -80
latmax = -80
lonmin = 360-48
lonmax = 360

loc_xmin = [-2e6,0.2e6,360-55,360-35,360-46]
loc_xmax = [2.5e6,1.5e6,360-45,360-30,360-35]
loc_ymin = [-2.5e6,-1.5e6,-74,-74.75,-78.1]
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

def TS_diagram(run,latS,lonW,startyr,endyr,z=0,zab=False,zall=True,plot_lines=True,
               seasonal=False,runcmp=False,savepath=savepath_nersc):

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
    nz = 60 
    
    T = np.zeros((nt,nz)) 
    S = np.zeros((nt,nz)) 
    if runcmp:
        T2 = np.zeros((nt,nz)) 
        S2 = np.zeros((nt,nz)) 
    
    placename = ( str(round(latS)) + 'S' + str(round(lonW)) + 'W' ) 
    zmid = zmidfrommesh(fmesh,cellidx=[idx])
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

    plt.colorbar(scalarMap,label=cbtitle)
    axTS.set_ylim([-2.1,0.2])
    axTS.set_xlim([33.5,34.8])
    axTS.set_ylabel(varlabel[vartitle.index('T')]) 
    axTS.set_xlabel(varlabel[vartitle.index('S')])

    lncolor = 'black'
    lw1 = 1 
    plt.plot([34.0, 34.0], [-1.85, 0.2],  ':', color=lncolor, linewidth=lw1)
    plt.plot([34.5, 34.5], [-1.86, -1.5], ':', color=lncolor, linewidth=lw1)
    plt.plot([34.0, 35.0], [-1.5, -1.5],  ':', color=lncolor, linewidth=lw1)
    plt.plot([34.0, 35.0], [0.0, 0.0],    ':', color=lncolor, linewidth=lw1)
    plt.plot([33.5, 35.0], [-0.0575*33.5+0.0901, -0.0575*35.0+0.0901], 
                                          ':', color=lncolor, linewidth=lw1)
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
                        zmid = zmidfrommesh(fmesh,cellidx=[idx])
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
        zmid = zmidfrommesh(fmesh,cellidx=[idx],vartype='scalar')
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

# PICK_TRANSECT
# option choose transect by 'coord' or 'contour'
#        if 'contour' need to provide latitude limits throgh latS and lonW
#        if 'coord' only need to provide contour endpoints
# contourvar  variable to choose contour from
# contourval  value of contour
# vartype 'scalar' or 'velocity
#
# Outputs:
#    cellidx
#    edgeidx
#    dist:  distance between each cell center and the cell center of the first cell
#    angle: mean orientation of transect from left to right
def pick_transect(option='coord',lat=[-76,-76],lon=[360-32,360-32],run='ISMF',
                  vartype='velocity',transect_name='',scope_name = 'frisEAcoast',
                  overwrite=False, savepath=savepath_nersc,append=False):
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])

    # import variables from file
    if vartype == 'scalar':
        latpt = fmesh.variables['latCell'][:]
        lonpt = fmesh.variables['lonCell'][:]
        xpt   = fmesh.variables['xCell'][:]
        ypt   = fmesh.variables['yCell'][:]
    elif vartype== 'velocity':
        latcell = fmesh.variables['latCell'][:]
        loncell = fmesh.variables['lonCell'][:]
        xcell = fmesh.variables['xCell'][:]
        ycell = fmesh.variables['yCell'][:]
        latpt = fmesh.variables['latEdge'][:]
        lonpt = fmesh.variables['lonEdge'][:]
        #idxpt = fmesh.variables['indexToEdgeID'][:]
        xpt   = fmesh.variables['xEdge'][:]
        ypt   = fmesh.variables['yEdge'][:]
   
    # the outcome of all these options is a list of x,y points with varying spacing and number 
    if option == 'coord':
        if transect_name == 'ryan_2017_trough':
            lat = [latmin,latmax]
            lon = [lonmin,lonmax]
        else:
            transect_name = (str(int(abs(lat[0]))) + 'S' + 
                         str(int(abs(lon[0]-360))) + 'W-' + 
                         str(int(abs(lat[1]))) + 'S' + 
                         str(int(abs(lon[1]-360))) + 'W' )
       #    dlat = 0.15 # at 30km resolution, distance between cells in deg
       #    dlon = 0.98
       #    # define plot location
       #    # define line of constant latitude
       #    if lat[0] == lat[1]:
       #        lat[0] = (lat[0] - dlat) * deg2rad
       #        lat[1] = (lat[1] + dlat) * deg2rad
       #    if lon[0] == lon[1]:
       #        lon[0] = (lon[0] - dlon) * deg2rad
       #        lon[1] = (lon[1] + dlon) * deg2rad
        #    # indices of transect
        #    logical_trans = ( (latcell> lat[0]) & 
        #                      (latcell< lat[1]) &
        #                      (loncell> lon[0]) & 
        #                      (loncell< lon[1])   )
        #    #logical_trans = ( (latpt> lat[0]*deg2rad) & 
        #    #                  (latpt< lat[1]*deg2rad) &
        #    #                  (lonpt> lon[0]*deg2rad) & 
        #    #                  (lonpt< lon[1]*deg2rad)   )
        #    idx = np.argmin(np.square(latcell[logical_trans]-lat[0])+
        #                    np.square(loncell[logical_trans]-lon[0]) )
        #    # sort by distance from the point that is furthest south
        #    # would be better to sort by distance from the point closest to latS[0],lonW[0]
        #    #idx1 = np.argmin(ypt[logical_trans])
        #    #dist = np.sqrt( np.square(ypt[logical_trans] - ypt[logical_trans][idx1]) + 
        #    #                np.square(xpt[logical_trans] - xpt[logical_trans][idx1])   )
        #    #idx = idx_trans[dist.argsort()]
        #    idx1 = np.argmin(ycell[idx])
        #    dist = np.sqrt( np.square(ycell[idx] - ycell[idx][idx1]) + 
        #                    np.square(xcell[idx] - xcell[idx][idx1])   )
        #    idx = idx[dist.argsort()]
        #    #print('number of cells:',str(len(idx)))
        #    print(latcell[idx])
        #    #ddist[0] = 0
        #    #ddist[1:] = dist[1:]-dist[:-1]
        #    #idx_max, = np.where(ddist>100e3)
        #    #if len(idx_max) > 0:
        #    #    print('large discontinuity found in transect')
        #    #    idx = idx[:idx_max[0]-1]# distance along transect
        #    #    print('new number of cells:',str(len(idx)))
        #elif select_cell == 'connecting':
        #    # NOT FUNCTIONAL

        #    #dlat = 0.3 # at 30km resolution, distance between cells in deg
        #    #dlon = 0.98
        #    #    lat[0] = (lat[0] - dlat) * deg2rad
        #    #    lat[1] = (lat[1] + dlat) * deg2rad
        #    #if lon[0] == lon[1]:
        #    #    lon[0] = (lon[0] - dlon) * deg2rad
        #    #    lon[1] = (lon[1] + dlon) * deg2rad

        #    # indices of transect
        #    #logical_trans = ( (latcell> (lat[0] - dlon) * deg2rad) & 
        #    #                  (latcell< (lat[0] + dlon) * deg2rad) &
        #    #                  (loncell> (lon[0] - dlon) * deg2rad) & 
        #    #                  (loncell< (lon[0] + dlon) * deg2rad) ) 
        #    #idx1_candidates = np.where(logical_trans)[0]
        #    #dist = np.sqrt( np.square(ycell[idx1_candidates] - 
        #    #                          ycell[idx][idx1]) + 
        #    #                np.square(xcell[idx1_candidates] - xcell[idx][idx1])   )
        #    dist = np.sqrt(np.square(latcell-lat[0]) + 
        #                   np.square(loncell-lon[0])) 
        #    idx1 = np.argmin(dist) 
        #    edgesOnCell = np.subtract(fmesh.variables['cellsOnCell'][idx,:],1)

        #edgesOnCell = np.subtract(fmesh.variables['edgesOnCell'][idx,:],1)
        
    elif option== 'zcontour':
        print('option ',option,' is not yet enabled') 
    
    elif option == 'by_index':
        datestr = '{0:04d}-{1:02d}'.format(98,1)
        #filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'
        #            .format(runpath[runname.index(run)]) 
        #            + datestr + '-01.nc')
        #f = netCDF4.Dataset(filename, 'r')
        
        if transect_name == 'trough_shelf':
            cellidx = np.subtract(cells_trough_shelf_lat,1)
        elif transect_name == 'trough_ice':
            cellidx = np.subtract(cells_trough_ice_lat,1)
        else:
            print('transect name not matched')
   
        #uCell = f.variables['timeMonthly_avg_velocityZonal'][0,idx,:]
        #vCell = f.variables['timeMonthly_avg_velocityMeridional'][0,idx,:]
    
    # TODO bit of code to identify cellidx for other cases    
    # constants
    # width of band used to draw points from, the larger the width, the more zig-zagging the line
    #res = 
    #if select_cell == 'within_halo':
    #    for i,_ in enumerate(x_transect):
    #        
    #    logical_trans = ( (xcell > x_transect[i]-res) & 
    #                      (xcell < lat[1]) &
    #                      (ycell > lon[0]) & 
    #                      (ycell < lon[1])   )

        x_transect = fmesh.variables['xCell'][cellidx]
        y_transect = fmesh.variables['yCell'][cellidx]
        edgesOnCell = np.subtract(fmesh.variables['edgesOnCell'][cellidx,:],1)
        verticesOnCell = np.subtract(fmesh.variables['verticesOnCell'][cellidx,:],1)
        verticesOnEdge = np.zeros((len(cellidx),7,2))
        for i in range(len(cellidx)):
            for j in range(7):
                verticesOnEdge[i,j,:] = (
                  fmesh.variables['verticesOnEdge'][edgesOnCell[i,j],:])
        
        x0 = xcell[cellidx[0]]
        y0 = ycell[cellidx[0]]
        x1 = xcell[cellidx[-1]]
        y1 = ycell[cellidx[-1]]
        m = (y1 - y0)/(x1 - x0)
        b = (x1*y0 - x0*y1)/(x1 - x0)
        angle = atan(1/m) #fmesh.variables['angleEdge'][edgesOnCell[i,j]])
        
        idx = [] 
        edgeidx = [] 
        dx = 5e3 
        for i,celli in enumerate(cellidx):
            for j,edge in enumerate(edgesOnCell[i,:]):
                ye = fmesh.variables['yEdge'][edge]
                xe = fmesh.variables['xEdge'][edge]
                ym = m*xe + b
                xm = (ye - b)/m
                #if ye > ym:
                if ( (xe < xm - dx) and (ye > y0) and (ye < y1) ):
                    #edge not in edgesOnCell[i+1,:]) ):
                    edgeidx.append(edgesOnCell[i,j])
                    idx.append(celli)
        cellidx = idx
        
        # TODO replace above method with that below
        # include all edges whose vertices have y-coordinates above or on a line 
        # connecting neighboring cell centers
#        for i,cell in enumerate(cellidx):
#            if cell == cellidx[-1]:
#                m = ( (y_transect[i] - y_transect[i-1]) /
#                      (x_transect[i] - x_transect[i-1]) )
#                b = ( (x_transect[i] * y_transect[i-1] - 
#                       x_transect[i-1] * y_transect[i]) /
#                      (x_transect[i] - x_transect[i-1]) )
#            else:
#                m = ( (y_transect[i+1] - y_transect[i]) /
#                      (x_transect[i+1] - x_transect[i]) )
#                b = ( (x_transect[i+1] * y_transect[i] - 
#                       x_transect[i] * y_transect[i+1]) /
#                      (x_transect[i+1] - x_transect[i]) )
#            for j,edge in enumerate(edgesOnCell[i,:]):
#                x2 = fmesh.variables['xVertex'][verticesOnEdge[i,j,:]
#                y2 = fmesh.variables['yVertex'][verticesOnEdge[i,j,:]
#                if ( ( y2[0] >= m*x2[0] + b) and 
#                     ( y2[1] >= m*x2[1] + b) and
#                     ( x2[1] >= m*x2[1] + b) and
#                     ( edge not in edgesOnCell[i+1,:]) ):
#                    edgeidx.append(edge)
                

        #    print('idxcell',str(i),':',str(ycell[i]),',',str(xcell[i]))
        #    shared_edges = []
        #    shared_verts = []
        #    xEdge = fmesh.variables['xEdge'][edgesOnCell[i,:]]))
        #    edge1 = np.argmin(xEdge)
        #       fmesh.variables[]verticesOnCell[i,:] = (
        #    vert1 = verticesOnCell[i,:]
        #    vert2 = verticesOnCell[i+1,:]
        #    edges1 = edgesOnCell[i,:]
        #    print('xedge:',str(fmesh.variables['xEdge'][edges1-1]))
        #    edges2 = fmesh.variables['edgesOnCell'][i+1,:]
        #    for j in vert1:
        #       if j in vert2:
        #          np.append(shared_verts,j)
        #    print('len(shared_edges) idx=',str(i),'=',str(len(shared_edges))) 
        #    for j in edges1:
        #       if j in edges2:
        #          np.append(shared_edges,j)
        #    print('len(shared_verts) idx=',str(i),'=',str(len(shared_verts))) 
        #for i in idx:
        #    print('idxcell',str(i),':',str(ycell[i]),',',str(xcell[i]))
        #    ii = fmesh.variables['indexToCellID'][:].index(i)
        #    for j,jj in enumerate(fmesh.variables['edgesOnCell'][ii,:]):
        #        print('edge ',str(j),':',str(fmesh.variables['yEdge'][jj]),',',str(fmesh.variables['xEdge'][jj]))
        #        plt.plot(fmesh.variables['yEdge'][jj],
        #                 fmesh.variables['xEdge'][jj],'r.',markersize=ms)
        #        k = fmesh.variables['cellsOnCell'][i,j]
        #        print('cell ',str(k),':',str(fmesh.variables['yCell'][k]),',',str(fmesh.variables['xEdge'][k]))
            #if (edges[1] in idx) and (edges[10] in idx):
            #   print('found duplicate')
            #   idx[idx.index(edges[10])] = nan
            #if (edges[6] in idx) and (edges[5] in idx):
            #   print('found duplicate')
            #   idx[idx.index(edges[5])] = nan
   
    #idx1 = np.argsort(fmesh.variables['yEdge'][edgeidx])
    #edgeidx = edgeidx[np.argsort(fmesh.variables['xEdge'][edgeidx])]
    dist = np.sqrt( np.square(fmesh.variables['yEdge'][edgeidx] - fmesh.variables['yEdge'][edgeidx[0]]) + 
                    np.square(fmesh.variables['xEdge'][edgeidx] - fmesh.variables['xEdge'][edgeidx[0]]) )
    
        # show profile line across cells
        #if not os.path.exists(savepath + 'bathy_' + placename + '.png'):
        #    fig1 = plt.figure()
        #    plt.plot(yCell[logical_N],     xCell[logical_N], 'k.')
        #    plt.plot(yCell[logical_trans], xCell[logical_trans], 'r.')
        #    plt.axis('equal')
        #    plt.savefig('grid_' + placename + '_' + datestr + '.png',dpi=set_dpi)
        #    plt.close()
        #    
        #    fig = plt.figure()
        #    cntr1 = plt.tricontourf(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
        #                            zmax[logical_N].flatten(), 20, cmap="viridis")
        #    plt.plot(yCell[logical_N],     xCell[logical_N], 'o', color = 'white', 
        #             markersize = 4, fillstyle = 'none')#, alpha = 0.5)
        #    cntr = plt.tricontour(yCell[logical_N].flatten(), xCell[logical_N].flatten(), 
        #                           zice[logical_N].flatten(), [-10], colors = 'k')
        #    plt.plot(yCell[idxsort_trans], xCell[idxsort_trans], 'k-')
        #    plt.axis('equal')
        #    cbar = plt.colorbar(cntr1)
        #    cbar.set_label('Depth (m)')    
        #    plt.savefig(savepath + 'bathy_' + placename + '.png',dpi=set_dpi)
        #    plt.close()
        
    # show profile line across cells
    if (not os.path.exists(savepath + 'bathy_' + transect_name + '.png')) or overwrite:
        # northern limit for subplots
        ms=1 
        #xcell = fmesh.variables['xCell'][:]
        #ycell = fmesh.variables['yCell'][:]
        #latcell = fmesh.variables['latCell'][:]
        scope_name1 = 'frisEAcoast'
        logical_N = ( (latcell < deg2rad*lat_N)                 & 
                      (xcell > loc_xmin[loc.index(scope_name1)]) & 
                      (xcell < loc_xmax[loc.index(scope_name1)]) & 
                      (ycell > loc_ymin[loc.index(scope_name1)]) & 
                      (ycell < loc_ymax[loc.index(scope_name1)]) )
        
        idx_scope = np.where(logical_N)[0]
        zmax_scope = np.multiply(-1.,fmesh.variables['bottomDepth'][idx_scope])
        icemask_scope  = fmesh.variables['landIceMask'][0,idx_scope]
        #zice_scope     = fmesh.variables['landIceDraft'][idx_scope]
        ycell_scope = ycell[idx_scope]
        xcell_scope = xcell[idx_scope]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cntr1 = plt.scatter(ycell_scope,xcell_scope,
                            s=loc_ptsize[loc.index(scope_name)], c=zmax_scope)
        #cntr1 = plt.tricontourf(ycell[idx_scope],xcell[idx_scope],zmax, 20, cmap="viridis")
        #plt.plot(ycell[idx_scope],xcell[idx_scope], 'o', color = 'white', 
        #         markersize = ms, fillstyle = 'none', alpha = 0.5)
        #cntr = plt.tricontour(ycell_scope,xcell_scope,
        #                      zice, [-10], colors = 'k',linewidth=lw1)
        plt.plot(ycell_scope[icemask_scope==1], xcell_scope[icemask_scope==1],
                 'o', color = 'white', 
                 markersize = 5*ms, fillstyle = 'none')
        plt.plot(ycell[cellidx], xcell[cellidx],
                 'o', color = 'black', 
                 markersize = ms, fillstyle = 'none')
        #plt.plot([y0,y1],[x0,x1], 'k-')

        cNorm  = Normalize(vmin=-1*pi, vmax=1*pi)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='jet')
        #for i in range(len(idx)):
        #    for j in range(6):
        #        colorVal = scalarMap.to_rgba(fmesh.variables['angleEdge'][edgesOnCell[i,j]])
        #        sc = plt.scatter(fmesh.variables['yEdge'][edgesOnCell[i,j]],
        #                    fmesh.variables['xEdge'][edgesOnCell[i,j]],s=ms/2,c=colorVal)
        for k in edgeidx:
            colorVal = scalarMap.to_rgba(fmesh.variables['angleEdge'][edgesOnCell[i,j]])
            sc = plt.plot([fmesh.variables['yVertex'][fmesh.variables['verticesOnEdge'][k,0]-1],
                           fmesh.variables['yVertex'][fmesh.variables['verticesOnEdge'][k,1]-1]],
                          [fmesh.variables['xVertex'][fmesh.variables['verticesOnEdge'][k,0]-1],
                           fmesh.variables['xVertex'][fmesh.variables['verticesOnEdge'][k,1]-1]],
                          'b-',linewidth=lw1)#marker='None',linestyle='-','k')
                #print(fmesh.variables['yVertex'][verticesOnEdge[i,j,1]])
                #print(fmesh.variables['xVertex'][verticesOnEdge[i,j,1]])
                #xv = [fmesh.variables['yVertex'][verticesOnEdge[i,j,0]],
                #      fmesh.variables['yVertex'][verticesOnEdge[i,j,1]]]
                #yv = [fmesh.variables['xVertex'][verticesOnEdge[i,j,0]],
                #      fmesh.variables['xVertex'][verticesOnEdge[i,j,1]]]
                #plt.plot(xv,yv,'b-',linewidth=0.5)
       #plt.plot(ypt[idx], xpt[idx], 'k.',markersize=ms)
        #    for i,j in enumerate(idx):
        #        plt.plot([yv1[i],yv2[i]],[xv1[i],xv2[i]], 'k-')
        #ax.set_xlabel('y (m)')
        #ax.set_ylabel('x (m)')
        plt.axis('equal')
        plt.xlim([loc_ymin[loc.index(scope_name)],loc_ymax[loc.index(scope_name)]])
        plt.ylim([loc_xmin[loc.index(scope_name)],loc_xmax[loc.index(scope_name)]])
        cbar = plt.colorbar(cntr1)
        #cbar = plt.colorbar(sc)
        cbar.set_label(r'Depth (m)')    
        print('save ','bathy_' + transect_name + '2')
        plt.savefig(savepath + 'bathy_' + transect_name + '2')
        plt.close()
    
    return cellidx, edgeidx, dist, angle


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
    zmid = zmidfrommesh(fmesh,cellidx=cellidx,vartype='scalar')
    zbottom = np.subtract(zmid,0.5*zh)
    ztop = np.add(zbottom,zh)
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
             lat=[latmin,latmax],lon=[lonmin,lonmax], 
             varlim = False, normal = False, zscale = 'linear', 
             run='ISMF', runcmp = False, runcmpname = 'ISMF-noEAIS',
             overwrite = False, ops = [''],
             save_transect_mean = False,
             savepath=savepath_nersc):
    
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
    depths   = fmesh.variables['refBottomDepth'][:]

    _,ymesh  = np.meshgrid(depths, dist)
    
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
    #zh       = fmesh.variables['layerThickness'][0,:]
    #icemask  = fmesh.variables['landIceMask'][idx]
    zice     = fmesh.variables['landIceDraft'][0,cellidx]
    bathymax = np.min(zmax) - 100
   
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
                                   str(varlim) +'.png') 

                # calculate z from depths
                zmesh = zmidfrommesh(fmesh, cellidx = cellidx, 
                                     vartype = 'scalar')
                
                if not overwrite and os.path.exists(image_filename):
                    print(image_filename + ' exists')
                    continue
                data     = f.variables[varname[vartitle.index(var)]][0,cellidx,:]
                #print('data_shape = ',np.shape(data))
                #print('ymesh_shape = ',np.shape(ymesh))
                #print('zmesh_shape = ',np.shape(zmesh))
                #print('data[0,:] = ',data[0,:])
                if runcmp: 
                    data2 = f2.variables[varname[vartitle.index(var)]][0,cellidx,:]
                
                if var == 'u' and normal:
                    u = data
                    v = f.variables[varname[vartitle.index('v')]][0,cellidx,:]
                    u_angle = np.arctan2(v,u)
                    transect_angle_normal = np.add(angle, pi/2)
                    u_norm = np.sqrt(np.add(np.square(u),np.square(v)))
                    data = np.multiply(u_norm, 
                                       np.cos(np.subtract(u_angle,
                                              transect_angle_normal)))
                # define plot location
                #if new:
                #    idx1 = np.argmin( np.square(yCell-lat_transect[0]) + 
                #                      np.square(xCell-lon_transect[0])   )
                #    idx2 = np.argmin( np.square(yCell-lat_transect[1]) + 
                #                      np.square(xCell-lon_transect[1])   )
                #    return
                #else:
                #    # define line of constant latitude
                #    if lat_transect[0] == lat_transect[1]:
                #        lat_transect[0] = lat_transect[0] - dlat
                #        lat_transect[1] = lat_transect[1] + dlat
                #    if lon_transect[0] == lon_transect[1]:
                #        lon_transect[0] = lon_transect[0] - dlon
                #        lon_transect[1] = lon_transect[1] + dlon
                #    # northern limit for subplots
                #    logical_N = (latCell < lat_N*deg2rad) & (xCell > 0)
                #    
                #    # indices of transect
                #    logical_trans = ( (latCell > lat_transect[0]*deg2rad) & 
                #                      (latCell < lat_transect[1]*deg2rad) &
                #                      (lonCell > lon_transect[0]*deg2rad) & 
                #                      (lonCell < lon_transect[1]*deg2rad)   )
                #    idx_trans = np.where(logical_trans)[0]
                #    
                #    idx1 = np.argmin(yCell[logical_trans])
                #    temp = np.sqrt( np.square(yCell[logical_trans] - yCell[logical_trans][idx1]) + 
                #                    np.square(xCell[logical_trans] - xCell[logical_trans][idx1])   )
                #    idxsort_trans = idx_trans[temp.argsort()]
                #    ysort_trans = yCell[idxsort_trans]
                #    xsort_trans = xCell[idxsort_trans]
                #    dist_trans = temp[temp.argsort()]
                #    dd_trans = np.zeros(dist_trans.shape)
                #    dd_trans[1:] = dist_trans[1:]-dist_trans[:-1]
                #    idx2, = np.where(dd_trans>100e3) 
                #    if len(idx2) > 0:
                #        dmax = dist_trans[idx2[0]-1]# distance along transect
                #    else:
                #        dmax = np.max(dist_trans)
               
                
                #data_trans = np.transpose(data[idxsort_trans,:])
                #data_trans_masked = np.transpose(data[idxsort_trans,:])
                #data_trans_zmasked = np.transpose(data[idxsort_trans,:])
                
                #zssh = zbottom[:,0] + zh[:,0]
                # Ideally, mask both bad data and ice shelf
                _,zicemesh = np.meshgrid(depths,zice)

                #data_masked = np.ma.masked_where( (data < bad_data) |
                #                                  (zmesh > zicemesh), 
                #                                  data)
                #data_trans_masked = np.ma.masked_where( (data_trans < bad_data) |
                #                                        (zmesh > zicemesh), data_trans)
                #data_masked = data
                data_masked = np.ma.masked_where( (data < bad_data) |
                        (data == bad_data2), data)
                for idx,i in enumerate(cellidx):
                    data_masked[idx,:] = np.ma.masked_where( 
                            (zmesh[idx,:] > zice[idx]), data[idx,:])
                    for jdx,j in enumerate(depths):
                        data_masked[idx,jdx] = np.ma.masked_where( 
                            (jdx > kmax[idx] ), data[idx,jdx])
                mask = np.ma.getmask(data_masked)
                #print(ymesh[~mask].flatten()) 
                #print(zmesh[~mask].flatten()) 
                #print(data_masked[~mask].flatten())
                #print(data_masked[0,:])
                #print(zmesh[0,:])
                #print(ssh[0])
                #print(zice[0])
                #print(zmax[0])
                #print(ymesh[0,:])
                
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
                if clevels[0] > np.min(data_masked[~mask].flatten()):
                    clevels = np.append(np.min(data_masked[~mask].flatten()),clevels)
                if clevels[-1] < np.max(data_masked[~mask].flatten()):
                    clevels = np.append(clevels,np.max(data_masked[~mask].flatten()))
                
                fig = plt.figure()
                
                if runcmp:
                    cmap1 = "cmo.balance"
                else:
                    cmap1 = varcmap[vartitle.index(var)]
                
                if runcmp:
                    data1_masked = data_masked
                    data_masked = np.subtract(data1_masked,data2[~mask])
                
                cntr2 = plt.tricontourf(np.transpose(ymesh[~mask].flatten()), 
                                        np.transpose(
                                           np.abs(zmesh[~mask].flatten())), 
                                        np.transpose(
                                           data_masked[~mask].flatten()), 
                                        levels=clevels,cmap=cmap1)
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
                    plt.ylim([1,np.abs(np.min(zmesh[~mask].flatten())-50)])
                ax.invert_yaxis()
                cbar = plt.colorbar()
                cbar.set_label(var)
                plt.xlabel('Distance (km)')
                plt.ylabel('Depth (m)')
                #plt.xlim([np.min(dist)/1e3,np.max(dist)/1e3])
                plt.title(runtitle[runname.index(run)] + ': ' + datestr)
                
                plt.savefig(image_filename,dpi=set_dpi)
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

def plot_stresscurl_t(var,run=['ISMF'],yrrange=[70,71],
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
    for i in run:
        headings.append(i+'_curl')
    wr.writerow(headings)
    for yr in range(yrrange[0],yrrange[1]):
        for mo in range(1,13):
            times = (yr+(mo-1.0)/12.0)
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            filein = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index(run)]) + datestr + '-01.nc'
            f = netCDF4.Dataset(filein, 'r')
            if runcmp:
               filein2 = '{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(runpath[runname.index('ISMF-noEAIS')]) + datestr + '-01.nc'
               f2 = netCDF4.Dataset(filein2, 'r')
            
            taux = f.variables['timeMonthly_avg_windStressZonal'][0,logical_N]
            tauy = f.variables['timeMonthly_avg_windStressMeridional'][0,logical_N]
            if runcmp:
               taux_cmp = f2.variables['timeMonthly_avg_windStressZonal'][0,logical_N]
               tauy_cmp = f2.variables['timeMonthly_avg_windStressMeridional'][0,logical_N]
            
            taux_i = interp.griddata((x,y),taux,(xi,yi),method='linear')
            tauy_i = interp.griddata((x,y),tauy,(xi,yi),method='linear')
            if runcmp:
                taux_i_cmp = interp.griddata((x,y),taux_cmp,(xi,yi),method='linear')
                tauy_i_cmp = interp.griddata((x,y),tauy_cmp,(xi,yi),method='linear')
            
            dtauy_dx = nan*np.ones((len(x_reg)-1,len(y_reg)-1))
            dtaux_dy = nan*np.ones((len(x_reg)-1,len(y_reg)-1))
            for i in range(1,len(x_reg)):
                dtauy_dx[i-1,:] = np.divide(tauy_i[i,1:] - tauy_i[i-1,1:],dy)
            for j in range(1,len(y_reg)):
                dtaux_dy[:,j-1] = np.divide(taux_i[1:,j] - taux_i[1:,j-1],dx)
            curl = np.nanmean(np.subtract(dtauy_dx,dtaux_dy))
            curl_t.append(curl)
            if runcmp:
                for i in range(1,len(x_reg)):
                    dtauy_dx[i-1,:] = np.divide(tauy_i_cmp[i,1:] - tauy_i_cmp[i-1,1:],dy)
                for j in range(1,len(y_reg)):
                    dtaux_dy[:,j-1] = np.divide(taux_i_cmp[1:,j] - taux_i_cmp[1:,j-1],dx)
                curl_cmp = np.nanmean(np.subtract(dtauy_dx,dtaux_dy))
                curl_t_cmp.append(curl_cmp)

            if runcmp:
               wr.writerow([yr,mo,times,curl,curl_cmp])
            else:
               wr.writerow([yr,mo,times,curl])

            if map_output:
                ms = 50
                fig = plt.figure()
                ax = fig.add_subplot(111)
                #plt.contourf(xi,yi,taux_i)
                plt.scatter(x,y,s=ms,c=tauy,edgecolors='k')
                plt.scatter(xi,yi,s=ms,c=tauy_i,edgecolors='r')
                plt.xlabel('xi',fontsize=16)
                plt.ylabel('yi',fontsize=16)
                plt.savefig(savepath+'tau_y_sc.png',dpi=100)
                plt.close(fig)
                fig = plt.figure()
                ax = fig.add_subplot(111)
                #plt.contourf(xi,yi,taux_i)
                plt.scatter(x,y,s=ms,c=taux,edgecolors='k')
                plt.scatter(xi,yi,s=ms,c=taux_i,edgecolors='r')
                plt.xlabel('xi',fontsize=16)
                plt.ylabel('yi',fontsize=16)
                plt.savefig(savepath+'tau_x_sc.png',dpi=100)
                plt.close(fig)
        
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.plot(times,curl_t,label=runname[0])
        #if runcmp:
        #   ax.plot(times,curl_t_cmp,label=runname[1])
        #ax.set_xlabel(r'Year',fontsize=fs)
        #ax.set_ylabel(r'Wind stress curl (N m$^{-3}$)',fontsize=fs)
        #ax.legend(loc=9,bbox_to_anchor=(0.15,-0.15),fontsize=fs)
        #plt.savefig(savepath+filename+'.png',dpi=100,bbox_inches='tight')
        #plt.close(fig)

    return

def plot_stresscurl_t_diff(filename,tlim=[9999.,9999.],
                      savepath=savepath_nersc):
    df = pandas.read_csv(savepath+filename+'.txt')
    t = df['decyear'][:]
    curl = df['ISMF_curl'][:]
    curl_cmp = df['ISMF-noEAIS_curl'][:]
    dcurl = np.subtract(curl,curl_cmp)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t,curl,'-b',label='ISMF')
    ax.plot(t,curl_cmp,'-k',label='ISMF-noEAIS')
    #ax.plot(t,dcurl)
    ax.set_xlabel(r'Year',fontsize=fs)
    ax.set_ylabel(r'Wind stress curl '+runname[0]+'-'+runname[1]+
                  '(N m$^{-3}$)',fontsize=fs)
    if tlim[0] != 9999.:
       ax.set_xlim(tlim)
       filename = filename + '_'  + str(int(tlim[0])) +'-'+ str(int(tlim[1]))
    ax.legend(loc=9,bbox_to_anchor=(0.15,-0.15),fontsize=fs)
    plt.savefig(savepath+filename+'_diff.png',dpi=100,bbox_inches='tight')
    plt.close(fig)
    #for i,row in data:
    #   for j,entry in enumerate(row):
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
    zmid = zmidfrommesh(fmesh,cellidx=idx)
    
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

#------------------------------------------------------------------------------
# ZMIDFROMMESH 
# -- solve for z values at the middle of the cell where prognostic variables are solved 
#
# Inputs:
#    fmesh  data structure for the mesh
#    zidx   (optional) indices for z level (list of integers)
#           if zidx=-1 (default), return the full z vector
#------------------------------------------------------------------------------
def zmidfrommesh(fmesh,zidx=[-1],cellidx=[],vartype='scalar'): 
    if vartype=='scalar':
        zh       = fmesh.variables['layerThickness'][0,cellidx]
        zmax     = fmesh.variables['bottomDepth'][cellidx]
    elif vartype== 'velocity':
        cell1 = np.subtract(fmesh.variables['cellsOnEdge'][cellidx,0],1)
        cell2 = np.subtract(fmesh.variables['cellsOnEdge'][cellidx,1],1)
        zh1   = fmesh.variables['layerThickness'][0,cell1,:]
        zh2   = fmesh.variables['layerThickness'][0,cell2,:]
        zh = (zh1 + zh2)/2
        zmax1 = fmesh.variables['bottomDepth'][cell1]
        zmax2   = fmesh.variables['bottomDepth'][cell2]
        zmax = (zmax1 + zmax2)/2
        #nz1 = fmesh.variables['maxLevelCell'][cell1]
        #nz2 = fmesh.variables['maxLevelCell'][cell2]
        #nz= nz1
        #for i in range(len(cell1)):
        #    nz[i] = np.min([nz1[i],nz2[i]])

    cells,zlevels = zh.shape
    # calculate z from depths
    #if zidx[0]<0:
    #    zidx = np.ones((cells,),dtype=bool)
    #zloc = zmax[zidx]
    zbottom  = np.zeros(zh.shape)
    for i in range(len(cellidx)):
        zbottom[i,-1] = np.multiply(-1,zmax[i])
        #print(zh[i,:])
        for j in range(zlevels-2,-1,-1):
            zbottom[i,j] = zbottom[i,j+1] + zh[i,j+1]
    zmid = zbottom + np.multiply(0.5,zh)
    ztop = zbottom + zh
    return zmid#,zbottom,ztop
