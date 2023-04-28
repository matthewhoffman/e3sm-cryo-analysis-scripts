#!/usr/bin/env python
'''
Script to plot melt time series from different runs
'''

import os
import glob
import pandas
import csv
from data_config import *
from extract_data import extract_tseries
from plot_config import *
from region_config import *
import pyproj
import xarray as xr
from pyremap import ProjectionGridDescriptor
import matplotlib.pyplot as plt

run_list = runname
year_range = [10, 100]
region = 'EAcoast'
#extract_tseries([run_list], ['land_ice_fw_flux'], year_range,
#                placename = 'fris', land_ice_mask = False)
var = 'land_ice_fw_flux'
compute_mpas_region_melt = False
compute_Rignot_melt = False
get_Rignot_melt_by_shelf = True
get_MPAS_melt_by_shelf = False
create_plot = True
ice_shelves_to_plot = ['Filchner', 'Ronne', 'EAshelves']
ice_shelves_titles = ['Filchner', 'Ronne', 'Eastern Weddell']
EAshelves = ['Brunt_Stancomb', 'Riiser-Larsen', 'Quar', 'Ekstrom']
#filename = f'{savepath}/bogus'
#if not os.path.exists(filename):
if compute_mpas_region_melt:
    extract_tseries(run_list, [var], year_range,
                    placename='EAcoast', land_ice_mask=True,
                    operation='area_mean')

if compute_Rignot_melt:
    obs_filename = f'{obspath}/Ocean/Melt/Rignot_2013_melt_rates_6000.0x6000.0km_10.0km_Antarctic_stereo.nc'
    obs_dataset = xr.open_dataset(obs_filename)
    lat = obs_dataset.lat.values
    lon = obs_dataset.lon.values
    meltRate = obs_dataset.meltRate.values
    idx_bool = ((lon < region_coordbounds[region_name.index(region),0,1] - 360)
                 & (lon > region_coordbounds[region_name.index(region),0,0] - 360))
    cellidx = np.asarray(idx_bool.nonzero(),dtype=int)[0,:]
    print(len(cellidx))
    # Each cell is 10km^2 so we do not need to weight by cell area
    meanMeltRate = np.nanmean(meltRate[cellidx])
    print(f'Computed mean melt: {meanMeltRate}')

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
    print(ice_shelf_names.keys())
    obsDict = {}  # dict for storing dict of obs data
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
            iceDict[shelfName] = {}
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
    # print(obsDict)
    meanMeltRate = np.zeros(len(obsFileNameDict.keys()))
    meanMeltRateUncertainty = np.zeros(len(obsFileNameDict.keys()))
    for i, obsName in enumerate(obsFileNameDict):
        print(obsName)
        cumArea = 0
        for shelfName in EAshelves:
            EADict = obsDict[obsName][shelfName]
            cumArea = cumArea + EADict['area']
        print(f'cumArea = {cumArea}')
        for shelfName in EAshelves:
            areaWeight = obsDict[obsName][shelfName]['area']/cumArea
            print(f'{shelfName}: melt={obsDict[obsName][shelfName]["meltRate"]}, area={obsDict[obsName][shelfName]["area"]}, weight={areaWeight}')
            meanMeltRate[i] = meanMeltRate[i] + areaWeight * obsDict[obsName][shelfName]['meltRate']
            meanMeltRateUncertainty[i] = meanMeltRateUncertainty[i] + areaWeight * obsDict[obsName][shelfName]['meltRateUncertainty']
        iceDict = obsDict[obsName]
        iceDict['EAshelves'] = {
            'area': cumArea,
            'meltFlux': 0,
            'meltFluxUncertainty': 0,
            'meltRate': meanMeltRate,
            'meltRateUncertainty': meanMeltRateUncertainty}
    obsDict[obsName] = iceDict
    print(f'Obs: {meanMeltRate[0]} +/- {meanMeltRateUncertainty[0]}')
    print(f'SS: {meanMeltRate[1]} +/- {meanMeltRateUncertainty[1]}')

if get_MPAS_melt_by_shelf:
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

if create_plot:
   # figure
   fig1 = plt.figure(1, facecolor='w', figsize=(4,8))
   nrow = len(ice_shelves_to_plot)
   ncol = 1
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
               
           axMelt.set_title(ice_shelves_titles[i], fontsize=10)
           axMelt.plot(times[mpas_melt != 0], mpas_melt[mpas_melt != 0], '-',
                       label=runtitle[runname.index(run)],
                       color=run_color[runname.index(run)],
                       linewidth=0.5
                      )
       boxHalfWidth = 10.
       obsName = 'Rignot et al. (2013)'
       obsShelf = obsDict[obsName][ice_shelf]
       #axMelt.errorbar(50,
       #             obsShelf['meltRate'],
       #             yerr=obsShelf['meltRateUncertainty'],
       #             fmt='*',
       #             color='orange', ecolor='orange',
       #             capsize=0,
       #             label=obsName)
       # plot a box around the error bar to make it more visible
       boxHalfHeight = obsShelf['meltRateUncertainty']
       print(boxHalfHeight)
       boxX = 50 + \
           boxHalfWidth * np.array([-1, 1, 1, -1, -1])
       boxY = obsShelf['meltRate'] + \
           boxHalfHeight * np.array([-1, -1, 1, 1, -1])
       axMelt.plot(boxX, boxY, '-', color='orange', linewidth=1,
                   label=obsName)
       obsName = 'Rignot et al. (2013) SS'
       obsShelf = obsDict[obsName][ice_shelf]
       #axMelt.errorbar(150,
       #             obsShelf['meltRate'],
       #             yerr=obsShelf['meltRateUncertainty'],
       #             fmt='*',
       #             color='dodgerblue', ecolor='dodgerblue',
       #             capsize=0,
       #             label=obsName)
       # plot a box around the error bar to make it more visible
       boxHalfHeight = obsShelf['meltRateUncertainty']
       print(boxHalfHeight)
       boxX = 150 + \
           boxHalfWidth * np.array([-1, 1, 1, -1, -1])
       boxY = obsShelf['meltRate'] + \
           boxHalfHeight * np.array([-1, -1, 1, 1, -1])

       axMelt.plot(boxX, boxY, '-', color='dodgerblue', linewidth=1,
                   label=obsName)
       axMelt.set_ylim([1e-1, 6])
       axMelt.set_yscale('log')
       if i+1 == len(ice_shelves_to_plot):
          axMelt.set_xlabel('Year', fontsize=10)
       else:
          axMelt.set_xticklabels('')
       if i == 0:
          h, l = axMelt.get_legend_handles_labels()
          axMelt.legend(h[0:6], l[0:6], frameon=False, fontsize=10, loc='upper right')
       axMelt.set_ylabel('Melt rate (m yr$^{-1}$)', fontsize=10)
       plt.savefig(f'{savepath}/mpas_melt_rates.png',dpi=set_dpi)
#print(melt_rate)
