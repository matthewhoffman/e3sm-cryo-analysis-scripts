#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

import os
import csv
import netCDF4
from pick_from_mesh import pick_from_region, pick_point
import numpy as np
import math
from extract_depths import zmidfrommesh
from plot_config import *
from data_config import *
from region_config import *

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
def extract_tseries(runlist,varlist,year_range,
                    placename = '', land_ice_mask=False,
                    lat=-9999,lon=-9999,
                    zrange=None, zeval=-9999, zab=False,
                    ztop_pyc = [False],zbottom_pyc = [False],
                    operation = 'mean',
                    overwrite=True, output_filename = '',
                    savepath=savepath): 

    if zab:
        m = 'mab'
    else:
        m = 'm'
    
    if len(ztop_pyc)<len(varlist):
        ztop_pyc = [False for i in varlist]
    if len(zbottom_pyc)<len(varlist):
        zbottom_pyc = [False for i in varlist]
    if output_filename == '':
        filename = ('_'.join(runlist) + '_' + 
                    ''.join(varlist) + '_' + placename)
        if zrange is not None:
            filename = filename + '_z{0:03d}-{1:03d}'.format(zrange[0], zrange[1]) + m
        filename = filename + '_t{0:03d}-{1:03d}'.format(year_range[0], year_range[1])
    else:
        filename = output_filename
    print('extract tseries',savepath + filename + '.txt')
    
    years = np.arange(year_range[0],year_range[1],1)
    months = np.arange(1,13,1)
    nt = len(years)*len(months)
    times = np.zeros((nt,))
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(runlist[0])])
    data = np.zeros((len(runlist),len(varlist),nt)) 
    if placename in region_is_point:
        lat = region_coordbounds[region_name.index(placename)][1,1]
        lon = region_coordbounds[region_name.index(placename)][0,0]
        idx,_ = pick_point(lat=lat,lon=lon,run=runlist[0],
                           plot_map=False,savepath=savepath)
        idx = [idx]
    elif lat != -9999:
        idx,_ = pick_point(lat=lat,lon=lon,run=runlist[0],
                           plot_map=False,savepath=savepath)
        idx = [idx]
    else:
        idx = pick_from_region(region=placename, run=runlist[0],
                               plot_map=True)

    #if 'unormal' in varlist:
    #    _,_,_,transect_angle = pick_transect(option='by_index',
    
#                                         run=run,transect_name = 'trough_shelf') 
    t=0
    colheadings = ['decyear']
    for j,run in enumerate(runlist):
        for i,var in enumerate(varlist):
            header = run+'_'+var
            if zeval != -9999:
                header = header + '_z' + str(int(zeval))
            if ztop_pyc[i]:
                header = header + '_abovepyc'
            if zbottom_pyc[i]:
                header = header + '_belowpyc'
            colheadings.append(header)
    
    if zrange is not None:
        kmax     = fmesh.variables['maxLevelCell'][idx]
        zmid,_,_ = zmidfrommesh(fmesh,cellidx=idx)
        zidx = np.zeros((len(idx),2),dtype=int)
        for i in range(0,len(idx)):
            if zrange[1] != -9999:
                if zab:
                    zeval = np.add(zmid[0][-1],zeval)
                zidx[i,:] = ([np.argmin(np.abs(np.subtract(zmid,zeval[0]))),
                         np.argmin(np.abs(np.subtract(zmid,zeval[1])))])
                if zidx[i,1] == zidx[i,0]:
                    zidx[i,1] += 1
                if zidx[i,1] < zidx[i,0]:
                    zidx[i,:] = [zidx[i,1],zidx[i,0]] 
            elif zeval != -9999:
                if zab:
                    zeval = np.add(zmid[0][-1],zeval)
                else:
                    zeval = -1.*zeval
                zidx[i,0] = np.argmin(np.abs(np.subtract(zmid[0],zeval)))
    
    for yr in years:
        print(yr)
        for mo in months:
            times[t] = yr+(mo-1.0)/12.0
       
            datestr = '{0:04d}-{1:02d}'.format(yr, mo)
            for j,run in enumerate(runlist):
                input_filename = ('{0}/mpaso.hist.am.timeSeriesStatsMonthly.'.format(
                            runpath[runname.index(run)]) 
                           + datestr + '-01.nc')
                if not os.path.exists(input_filename):
                   print('Output file for {} does not exist'.format(run))
                   print(f'looking at {runpath[runname.index(run)]}')
                   print(f'for file {input_filename}')
                   data[j,:,t] = math.nan
                   continue
                f = netCDF4.Dataset(input_filename, 'r')
                if any(ztop_pyc) or any(zbottom_pyc):
                    z_idx = np.zeros((len(idx)),dtype=int)
                    zcol_mean = np.zeros((len(idx)))# make idx a list
                    T = f.variables[varname[vartitle.index('T')]][0,idx,:]
                    S = f.variables[varname[vartitle.index('S')]][0,idx,:]
                    h = fmesh.variables['layerThickness'][0,idx,:]
                    for idx_i,_ in enumerate(idx):
                        zpyc = z_pycnocline(zmid[idx_i,:kmax[idx_i]],
                                            T   [idx_i,:kmax[idx_i]],
                                            S   [idx_i,:kmax[idx_i]],
                                            zmin = -500.,cellidx=idx[idx_i])
                        z_idx[idx_i] = int(np.argmin(np.abs(np.subtract(zmid[idx_i,:kmax[idx_i]],zpyc))))
                for i,var in enumerate(varlist):
                    if var in surfvar:
                        if operation == 'mean':
                           data[j,i,t] = np.mean(f.variables[varname[vartitle.index(var)]][0,idx])
                           print(f'{yr}-{mo} mean var: {data[j,i,t]}')
                        elif operation == 'area_sum':
                           data[j,i,t] = np.sum(np.multiply(f.variables[varname[vartitle.index(var)]][0,idx],
                                                            fmesh.variables['areaCell'][idx]))
                        elif operation == 'area_mean':
                           data[j,i,t] = np.sum(np.multiply(f.variables[varname[vartitle.index(var)]][0,idx],
                                                            fmesh.variables['areaCell'][idx]))
                           A = fmesh.variables['areaCell'][idx]
                           if land_ice_mask:
                               mask = f.variables[varname[vartitle.index(var)]][0,idx] != 0
                               A = A[mask]
                               data[j,i,t] = data[j,i,t]/np.sum(A)
                        else:
                           data[j,i,t] = f.variables[varname[vartitle.index(var)]][0,idx]
                    else:
                        if ztop_pyc[i] or zbottom_pyc[i]:
                            var = f.variables[varname[vartitle.index(var)]][0,idx,:]
                            for idx_i,_ in enumerate(idx):
                                if np.isnan(z_idx[idx_i]) or z_idx[idx_i] == 0 or z_idx[idx_i] == kmax[idx_i]:
                                    zcol_mean[idx_i] = -9999
                                else:
                                    if ztop_pyc[i]:
                                        idx_top = 0
                                        #idx_top = np.minimum(z_idx-1,
                                        #            np.argmin(np.abs(np.subtract(zmid[idx_i,:kmax[idx_i]],zrange[0]))))
                                        zcol_mean[idx_i] = np.divide(
                                                             np.sum(np.multiply(
                                                               h[idx_i,idx_top:z_idx[idx_i]],
                                                               var[idx_i,idx_top:z_idx[idx_i]])),
                                                               np.sum(h[idx_i,idx_top:z_idx[idx_i]]))
                                    elif zbottom_pyc[i]:
                                        idx_bottom = kmax[idx_i]
                                        #idx_bottom = np.minimum(kmax[idx_i],
                                        #               np.argmin(np.abs(np.subtract(zmid[idx_i,:kmax[idx_i]],zrange[1]))))
                                        zcol_mean[idx_i] = np.divide(
                                                             np.sum(np.multiply(
                                                               h       [idx_i,z_idx[idx_i]:idx_bottom],
                                                               var     [idx_i,z_idx[idx_i]:idx_bottom])),
                                                               np.sum(h[idx_i,z_idx[idx_i]:idx_bottom]))
                            data[j,i,t] = np.nanmean(zcol_mean[zcol_mean != -9999])
                        elif zrange[1] != -9999:
                            data[j,i,t] = np.mean(f.variables[varname[vartitle.index(var)]]
                                                             [0,idx,zidx[:,0]:zidx[:,1]]       )
                        elif zeval != -9999:
                            data[j,i,t] = f.variables[varname[vartitle.index(var)]][0,idx,zidx[:,0]]
                f.close()

            t += 1
    
    if overwrite:
        flag='w+'
    else:
        flag='a+'
    table_file = open(f'{savepath}/{filename}.txt',flag)
    wr = csv.writer(table_file,dialect='excel')
    wr.writerow(colheadings)
    rowentries = np.zeros((len(varlist)*len(runlist)+1))
    for i,t in enumerate(times):
        rowentries[0] = t
        rowentries[1:] = data[:,:,i].flatten()
        wr.writerow(rowentries)

    return 

def extract_MPAS_melt_by_shelf(run_list, ice_shelves_to_plot):
    for run in run_list:
        filelist = glob.glob(f'{analysispath[runname.index(run)]}/timeseries/iceShelfFluxes/iceShelfFluxes*.nc')
        decYearMpas = np.zeros((len(filelist*12)))
        meltRateMpas = np.zeros((len(ice_shelves_to_plot), len(filelist*12)))
        f = xr.open_dataset(filelist[0])
        regionNames = f.regionNames.values.tolist()
        k = 0
        for j, infile in enumerate(sorted(filelist)):
            f = xr.open_dataset(infile)
            Year = int(infile[-7:-3])
            print(Year)
            for i, shelfName in enumerate(ice_shelves_to_plot[:-1]):
                data = f['meltRates'].isel(nRegions=regionNames.index(shelfName))
                #meltRateMpas[i, j] = data.mean(dim='Time') # each time index holds 1 year
                if np.shape(data.values)[0] > 12:
                    print(f'Skipping {infile}')
                    continue
                meltRateMpas[i, k:k+12] = data.values # each time index holds 1 year
            decYearMpas[k:k+12] = [Year + i/12 for i in range(12)]
            k = k+12
        ds = xr.Dataset()
        ds['meltRate'] = xr.DataArray(meltRateMpas, dims=('iceShelf', 'time'))
        ds['times'] = xr.DataArray(decYearMpas, dims=('time'))
        ds.to_netcdf(f'{savepath}/{run}_Filchner-Ronne-MPAS-melt-rates.nc')
    return
