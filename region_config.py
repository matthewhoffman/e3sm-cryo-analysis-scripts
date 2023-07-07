#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

import gsw
import numpy as np

cells_trough_ice_lat = [33000,77033,11788,129098,49348] # index is value - 1, located just under ice
cells_trough_shelf_lat = [165785,166563,130260,191987,85569] # index is value - 1, located just under ice
#cells_trough_shelf_lat = [224105,165785,166563,130260,191987] # index is value - 1, located just under ice
cells_trough_shelf2_lat = [152383,144943,67352] # index is value - 1, located just under ice

latmin = -80
latmax = -80
lonmin = 360-48
lonmax = 360

region_name = ['frisEAcoast','fris','gyre_interior',
               'trough_shelf',
               'etrough_crossshelf','trough_crossshelf','wtrough_crossshelf',
               'trough_ice','M31W_W','M31W','N31W','S4E','wedwang','wedwang_c',
               'wed_pyc_west','wed_pyc_filchner_trough','wed_pyc_brunt',
               'wed_pyc_Ryan','wed_pyc_Ryan_shallow','wed_pyc_Ryan_shelf',
               'EAcoast']
# M31W and S4E from Ryan et al. 2017 JGR Oceans
region_is_point = ['S4E','M31W','M31W_W','N31W']

wedwang_c_btr_flux_obs = [-6,-2]
gyre_interior_drho_obs = [0.2,0.4]
wed_pyc_Ryan_dzpyc_obs = [450]
M31W_T_obs = [-2,-1]
M31W_S_obs = [34.4,34.65]
M31W_rho_obs = [gsw.rho(34.5,-1.1,10),gsw.rho(34.6,-2,10)]

region_title = ['' for i in region_name] # loc_name
#TODO fix allocation
region_coordbounds = np.zeros((len(region_name),2,2)) 
region_xybounds    = np.zeros((len(region_name),2,2)) 
region_zbounds     = np.zeros((len(region_name),2)) 
for i in region_name:
    region_zbounds[region_name.index(i)] = [-9999.,20.]

region_title[region_name.index('EAcoast')] = 'Eastern Weddell Ice Shelves'
region_title[region_name.index('wedwang_c')] = 'Weddell Gyre'
region_title[region_name.index('wed_pyc_west')] = 'Eastern Weddell Shelf'
region_title[region_name.index('wed_pyc_Ryan_shelf')] = 'on-shelf'
region_title[region_name.index('wed_pyc_Ryan_shallow')] = 'on-shelf'
region_title[region_name.index('gyre_interior')] = 'off-shelf'
region_title[region_name.index('trough_shelf')] = 'trough, shelf break'
region_title[region_name.index('trough_shelf')] = 'trough, shelf break'
region_title[region_name.index('trough_ice')] = 'trough, under ice shelf front'
region_title[region_name.index('M31W_W')] = 'trough, west'
region_title[region_name.index('M31W')] = 'M31W'
region_title[region_name.index('N31W')] = 'WDW, shelf break'
region_title[region_name.index('S4E')] = 'S4E'

region_coordbounds[region_name.index('M31W_W')] = ([360-37,360-37],
                                                 [-76.046,-76.046])
region_coordbounds[region_name.index('M31W')] = ([360-30.994,360-30.994],
                                                 [-76.046,-76.046])
region_coordbounds[region_name.index('N31W')] = ([360-30.994,360-30.994],
                                                 [-72.8,-72.8])
region_coordbounds[region_name.index('S4E')] = ([360-30.58,360-30.58],
                                                [-74.62,-74.62])
region_xybounds   [region_name.index('EAcoast')]      = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('EAcoast')]      = [332,350],[-180,0]
#region_zbounds    [region_name.index('EAcoast')]      = [-5000,10]
region_xybounds   [region_name.index('frisEAcoast')]  = [0e6,2e6],[-2e6,0] 
region_coordbounds[region_name.index('frisEAcoast')]  = [300,360],[-180,-60]
region_xybounds   [region_name.index('fris')]         = [0.2e6,1.5e6],[-1.5e6,-0.4e6] 
region_coordbounds[region_name.index('fris')]         = [-1,361],[-181,181]
region_xybounds   [region_name.index('wedwang')]      = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wedwang')]      = [305,315],[-74,-65] 
region_xybounds   [region_name.index('wedwang_c')]    = [1.25e6,1.75e6],[-0.75e6,-0.75e6] 
region_coordbounds[region_name.index('wedwang_c')]    = [330,335],[-75,-72]               
region_zbounds    [region_name.index('wedwang_c')]    = [-3250,-750]
region_xybounds   [region_name.index('wed_pyc_west')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wed_pyc_west')] = [300,310],[-74,-72] 
region_xybounds   [region_name.index('wed_pyc_brunt')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wed_pyc_brunt')] = [330,340],[-75,-71] 
region_xybounds   [region_name.index('wed_pyc_filchner_trough')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wed_pyc_filchner_trough')] = [320,330],[-75,-71] 
region_xybounds   [region_name.index('wed_pyc_Ryan')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wed_pyc_Ryan')] = [330,360],[-75,-69] 
region_xybounds   [region_name.index('trough_shelf')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('trough_shelf')] = [323,332],[-74.75,-74.75]
region_xybounds   [region_name.index('trough_ice')]   = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('trough_ice')]   = [314,325],[-78.1,-78.1]
region_xybounds   [region_name.index('gyre_interior')]= region_xybounds[region_name.index('wed_pyc_Ryan')]
region_coordbounds[region_name.index('gyre_interior')]= [330,360],[-75,-68]
region_zbounds    [region_name.index('gyre_interior')]= [-4500,-3500]
#region_zbounds    [region_name.index('gyre_interior')]= [-8000,-4500]
region_xybounds   [region_name.index('wtrough_crossshelf')] = [1.25e6,1.75e6],[-0.75e6,-0.75e6] 
region_coordbounds[region_name.index('wtrough_crossshelf')] = [320,329],[-76,-72]
region_xybounds   [region_name.index('trough_crossshelf')] = [1.25e6,1.75e6],[-0.75e6,-0.75e6] 
region_coordbounds[region_name.index('trough_crossshelf')] = [327,335],[-76,-72]
region_xybounds   [region_name.index('etrough_crossshelf')] = [1.25e6,1.75e6],[-0.75e6,-0.75e6] 
region_coordbounds[region_name.index('etrough_crossshelf')] = [332,338],[-76,-72]
region_xybounds   [region_name.index('wed_pyc_Ryan_shallow')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wed_pyc_Ryan_shallow')] = [330,360],[-75,-69] 
region_zbounds    [region_name.index('wed_pyc_Ryan_shallow')]= [-500,0]
region_xybounds   [region_name.index('wed_pyc_Ryan_shelf')] = [-6.5e6,6.5e6],[-6.5e6,6.5e6] 
region_coordbounds[region_name.index('wed_pyc_Ryan_shelf')] = [330,360],[-75,-69] 
region_zbounds    [region_name.index('wed_pyc_Ryan_shelf')]= [-2000,-500]

loc_ptsize = [40 for i in region_name]

