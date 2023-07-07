#!/usr/bin/env python
'''
Script to plot melt time series from different runs
'''

import os
import glob
import pandas
import csv
from data_config import *
from extract_data import extract_tseries, extract_MPAS_melt_by_shelf
from plot_config import *
from region_config import *
import pyproj
import xarray as xr
from pyremap import ProjectionGridDescriptor
import matplotlib.pyplot as plt

run_list = runname
year_range = [10, 100]
region = 'EAcoast'
var = 'land_ice_fw_flux'
compute_mpas_region_melt = False
compute_Rignot_melt = False
get_Adusumilli_melt_by_shelf = True
get_Rignot_melt_by_shelf = True
get_MPAS_melt_by_shelf = False
create_plot = True
ice_shelves_to_plot = ['Filchner', 'Ronne', 'EAshelves']
ice_shelves_titles = ['Filchner', 'Ronne', 'Eastern Weddell']
EAshelves = ['Brunt_Stancomb', 'Riiser-Larsen', 'Quar', 'Ekstrom']
obsDict = {}  # dict for storing dict of obs data
iceDict = {}  # dict for storing dict of shelf data

if compute_mpas_region_melt:
    extract_tseries(run_list, [var], year_range,
                    placename='EAcoast', land_ice_mask=True,
                    operation='area_mean')

if get_Adusumilli_melt_by_shelf:
    obsName = 'Adusumilli et al. (2020)'
    obsDict[obsName] = {}
    ds = xr.open_dataset(f'{obspath}/Ocean/Melt/Adusumilli/Adusumilli_2020_iceshelf_melt_rates_2010-2018_v0.20230504.iceShelves20200621.nc')
    region_names = ds.regionNames
    mean_melt_rate = ds.meanMeltRate.values
    region_area = ds.meltArea.values
    melt_rate_sd = ds.meltRateUncertainty.values
    print(np.shape(mean_melt_rate))
    for region in ice_shelves_to_plot:
        if region == 'EAshelves':
            mean_melt_region = 0.
            melt_sdsq_area_weighted = 0.
            cum_area = 0.
            cum_area_sq = 0.
            for ice_shelf in EAshelves:
                idx = np.where(region_names==ice_shelf)
                cum_area += region_area[idx]
                cum_area_sq += region_area[idx]**2.
            for ice_shelf in EAshelves:
                idx = np.where(region_names==ice_shelf)
                area_weight = region_area[idx] / cum_area
                area_weight_sq = region_area[idx]**2. / cum_area_sq
                mean_melt_region += mean_melt_rate[idx] * area_weight
                melt_sdsq_area_weighted += melt_rate_sd[idx]**2. * area_weight_sq
                #print(f'{region_area[idx]},{cum_area}')
                #print(f'{area_weight},{region_area[idx]**2. / cum_area_sq}')
                #print(f'{ice_shelf}={mean_melt_rate[idx]}+-{melt_rate_sd[idx]}')
            melt_sd_region = melt_sdsq_area_weighted**0.5
        else:
            idx = np.where(region_names==region)
            mean_melt_region = mean_melt_rate[idx]
            melt_sd_region = melt_rate_sd[idx]
        print(f'{region}={mean_melt_region}+-{melt_sd_region}')
        iceDict[region] = {
            'area': region_area[idx],
            'meltFlux': 0.0,
            'meltFluxUncertainty': 0.0,
            'meltRate': mean_melt_region,
            'meltRateUncertainty': melt_sd_region}
        obsDict[obsName] = iceDict
        
if get_Rignot_melt_by_shelf:
    Rignot_data = 'Ocean/Melt/Rignot_2013_melt_rates.csv'
    # process obs
    obsFileNameDict = {'Rignot et al. (2013)':
                       'Rignot_2013_melt_rates_20201117.csv',
                       'Rignot et al. (2013) SS':
                       'Rignot_2013_melt_rates_SS_20201117.csv'}
    
    
    # former is Rignot, latter is MPAS
    ice_shelf_names = {'Filchner':'Filchner',
                       'Ronne':'Ronne',
                       'Brunt_Stancomb':'Brunt_Stancomb',
                       'Stancomb-Wills':None,
                       'Riiser-Larsen':None,
                       'Quar':None,
                       'Ekstrom':None
                      }
    # for ice_shelf_to_plot in ice_shelves_to_plot:
    for obsName in obsFileNameDict:
        obsFileName = f'{obspath}/Ocean/Melt/{obsFileNameDict[obsName]}'
        obsDict[obsName] = {}
        obsFile = csv.reader(open(obsFileName, 'rU'))
        iceDict = {}  # dict for storing dict of shelf data
        next(obsFile, None)  # skip the header line
        for line in obsFile:  # some later useful values commented out
            shelfName = line[0]
            if shelfName not in ice_shelf_names.keys():
                continue
            meltFlux = float(line[2])
            meltFluxUncertainty = float(line[3])
            meltRate = float(line[4])
            meltRateUncertainty = float(line[5])
            area = meltFlux * 1e9 / (meltRate * rho_ice * 907.185) 
            # build dict of obs. keyed to filename description
            # (which will be used for plotting)
            iceDict[shelfName] = {
                'area': area,
                'meltFlux': meltFlux,
                'meltFluxUncertainty': meltFluxUncertainty,
                'meltRate': meltRate,
                'meltRateUncertainty': meltRateUncertainty}
        obsDict[obsName] = iceDict

    meanMeltRate = np.zeros(len(obsFileNameDict.keys()))
    meanMeltRateUncertainty = np.zeros(len(obsFileNameDict.keys()))
    for i, obsName in enumerate(obsFileNameDict):
        print(obsName)
        cumArea = 0.
        cumAreaSq = 0.
        for shelfName in EAshelves:
            EADict = obsDict[obsName][shelfName]
            cumArea += EADict['area']
            cumAreaSq += EADict['area']**2.
        for shelfName in EAshelves:
            EADict = obsDict[obsName][shelfName]
            areaWeight = EADict['area']/cumArea
            areaWeightSq = EADict['area']**2. / cumAreaSq
            print(f'{shelfName}: melt={EADict["meltRate"]}, area={EADict["area"]}, weight={areaWeight}')
            meanMeltRate[i] += areaWeight * obsDict[obsName][shelfName]['meltRate']
            meanMeltRateUncertainty[i] += areaWeightSq * obsDict[obsName][shelfName]['meltRateUncertainty']**2.
            print(f'{shelfName}: melt_sd={EADict["meltRateUncertainty"]**2}, weight={areaWeightSq}')
        iceDict = obsDict[obsName]
        iceDict['EAshelves'] = {
            'area': cumArea,
            'meltFlux': 0,
            'meltFluxUncertainty': 0,
            'meltRate': meanMeltRate[i],
            'meltRateUncertainty': meanMeltRateUncertainty[i]**0.5}
    obsDict[obsName] = iceDict

if get_MPAS_melt_by_shelf:
    extract_MPAS_melt_by_shelf(run_list, ice_shelves_to_plot)

if create_plot:
    # figure
    fig1 = plt.figure(1, facecolor='w', figsize=(6,8))
    nrow = len(ice_shelves_to_plot)
    ncol = 1
    axis_label = ['a', 'b', 'c']
    for i, ice_shelf in enumerate(ice_shelves_to_plot):
        axMelt = fig1.add_subplot(nrow, ncol, i+1)
        for run in run_list:
            if ice_shelf == 'Filchner':
                mpas_filename = f'{savepath}/{run}_Filchner-Ronne-MPAS-melt-rates.nc'
                mpas_dataset = xr.open_dataset(mpas_filename)
                mpas_melt = mpas_dataset.meltRate.values[0,:]
                times = mpas_dataset.times.values
            elif ice_shelf == 'Ronne':
                mpas_filename = f'{savepath}/{run}_Filchner-Ronne-MPAS-melt-rates.nc'
                mpas_dataset = xr.open_dataset(mpas_filename)
                mpas_melt = mpas_dataset.meltRate.values[1,:]
                times = mpas_dataset.times.values
            elif ice_shelf == 'EAshelves':
                mpas_filename = f'{savepath}/CGM-UIB_CGM-DIB_VGM-DIB_land_ice_fw_flux_EAcoast_t001-200.txt'
                df = pandas.read_csv(mpas_filename)
                times = df['decyear'][:]
                melt_flux = df[f'{run}_{var}'][:]
                # the output is now in units of kg s^-1 m^-2
                mpas_melt = melt_flux * s_to_yr / rho_ice
                
            axMelt.set_title(f'{axis_label[i]}. {ice_shelves_titles[i]}', loc='left')
            axMelt.plot(times[mpas_melt != 0], mpas_melt[mpas_melt != 0], '-',
                        label=runtitle[runname.index(run)],
                        color=run_color[runname.index(run)],
                        linewidth=0.5
                       )
        boxHalfWidth = 10.
        obsName = 'Adusumilli et al. (2020)'
        obsShelf = obsDict[obsName][ice_shelf]
        # plot a box around the error bar to make it more visible
        boxHalfHeight = obsShelf['meltRateUncertainty']
        boxX = 50 + \
            boxHalfWidth * np.array([-1, 1, 1, -1, -1])
        boxY = obsShelf['meltRate'] + \
            boxHalfHeight * np.array([-1, -1, 1, 1, -1])
        axMelt.plot(boxX, boxY, '-', color='red', linewidth=1,
                    label=obsName)
        obsName = 'Rignot et al. (2013)'
        obsShelf = obsDict[obsName][ice_shelf]
        # plot a box around the error bar to make it more visible
        boxHalfHeight = obsShelf['meltRateUncertainty']
        boxX = 100 + \
            boxHalfWidth * np.array([-1, 1, 1, -1, -1])
        boxY = obsShelf['meltRate'] + \
            boxHalfHeight * np.array([-1, -1, 1, 1, -1])
        axMelt.plot(boxX, boxY, '-', color='orange', linewidth=1,
                    label=obsName)
        obsName = 'Rignot et al. (2013) SS'
        obsShelf = obsDict[obsName][ice_shelf]
        # plot a box around the error bar to make it more visible
        boxHalfHeight = obsShelf['meltRateUncertainty']
        boxX = 150 + \
            boxHalfWidth * np.array([-1, 1, 1, -1, -1])
        boxY = obsShelf['meltRate'] + \
            boxHalfHeight * np.array([-1, -1, 1, 1, -1])
        axMelt.plot(boxX, boxY, '-', color='dodgerblue', linewidth=1,
                    label=obsName)
        axMelt.set_xlim([0,200])
        #axMelt.set_ylim([-3e-2, 6])
        if i == 0:
           h, l = axMelt.get_legend_handles_labels()
           axMelt.legend(h[0:6], l[0:6], frameon=False, fontsize=10,
                         loc='upper right')
        if i+1 == len(ice_shelves_to_plot):
           axMelt.set_xlabel('Year')#, fontsize=10)
        else:
           axMelt.set_xticklabels('')
        axMelt.set_ylabel('Melt rate (m yr$^{-1}$)')
        plt.savefig(f'{savepath}/ice-shelf-melt-timeseries.png',dpi=set_dpi)
