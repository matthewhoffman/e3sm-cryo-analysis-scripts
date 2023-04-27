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
from math import pi,sqrt,nan,atan,cos,floor,ceil
import pandas
from extract_depths import zmidfrommesh
from plot_config import *

def pick_point(lat=-9999,lon=-9999,
               run='ISMF',placename = '',
               vartype='velocity',transect_name='',plot_region = 'frisEAcoast',
               plot_map=False, overwrite=False, savepath=savepath_nersc):
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])
    
    if lat == -9999:
        lat = region_coordbounds[region_name.index(placename)][1,1]
        lon = region_coordbounds[region_name.index(placename)][0,0]
    
    # import variables from file
    latCell = fmesh.variables['latCell'][:]
    lonCell = fmesh.variables['lonCell'][:]
    xCell   = fmesh.variables['xCell'][:]
    yCell   = fmesh.variables['yCell'][:]
    if vartype== 'velocity':
        latpt = fmesh.variables['latEdge'][:]
        lonpt = fmesh.variables['lonEdge'][:]
        #idxpt = fmesh.variables['indexToEdgeID'][:]
        xpt   = fmesh.variables['xEdge'][:]
        ypt   = fmesh.variables['yEdge'][:]
   
    # the outcome of all these options is a list of x,y points with varying spacing and number 
    geod = pyproj.Geod(ellps='WGS84')
    dlat = 0.3 # at 30km resolution, distance between cells in deg
    dlon = 0.98
    
    if placename == '':
        location_name = (str(int(abs(lat))) + 'S' + 
                     str(int(abs(lon-360))) + 'W')
    else:
        location_name = placename
 
    candidate_bool = ( (latCell > (lat-dlat)*deg2rad )    & 
                       (latCell < (lat+dlat)*deg2rad )    &
                       (lonCell > (lon-dlon)*deg2rad )    & 
                       (lonCell < (lon+dlon)*deg2rad )  )
    candidate_idx = np.asarray(candidate_bool.nonzero(),dtype=int)[0,:]
    _,_,distance_from_point = geod.inv(
                              [lon for j in candidate_idx],
                              [lat for j in candidate_idx],
                              lonCell[candidate_idx]/deg2rad,
                              latCell[candidate_idx]/deg2rad)
    cellidx = int(candidate_idx[np.argmin(distance_from_point)])

    if plot_map:
        fig = plt.figure()
        idx = pick_from_region(region=plot_region,run=run)
        zmax     = np.multiply(-1,fmesh.variables['bottomDepth'][idx])
        zice     = fmesh.variables['landIceDraft'][0,idx]
        
        cntr1 = plt.tricontourf(yCell[idx].flatten(), xCell[idx].flatten(),
                                zmax.flatten(), 20, cmap="viridis")
        #plt.plot(yCell[idx],     xCell[idx], 'o', color = 'white',
        #         markersize = 4, fillstyle = 'none')#, alpha = 0.5)
        plt.plot(yCell[cellidx], xCell[cellidx], 'o', color = 'red',
                 markersize = 4)
        #cntr = plt.tricontour(yCell[idx].flatten(), xCell[idx].flatten(),
        #                       zice.flatten(), [-10], colors = 'k')
        #if varlim:
        #    plt.clim([varmin[vartitle.index('z')], varmax[vartitle.index('z')]])
        plt.axis('equal')
        cbar = plt.colorbar(cntr1)
        cbar.set_label('Depth (m)')
        plot_filename = 'bathy_' + location_name
        print(plot_filename)
        plt.savefig(savepath + plot_filename + '.png',dpi=set_dpi)
        plt.close()

    return cellidx,location_name

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
    latcell = fmesh.variables['latCell'][:]
    loncell = fmesh.variables['lonCell'][:]
    xcell   = fmesh.variables['xCell'][:]
    ycell   = fmesh.variables['yCell'][:]
    zbottom = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    #if vartype == 'scalar':
        #latpt = fmesh.variables['latCell'][:]
        #lonpt = fmesh.variables['lonCell'][:]
    if vartype== 'velocity':
        latpt = fmesh.variables['latEdge'][:]
        lonpt = fmesh.variables['lonEdge'][:]
        #idxpt = fmesh.variables['indexToEdgeID'][:]
        xpt   = fmesh.variables['xEdge'][:]
        ypt   = fmesh.variables['yEdge'][:]
   
    # the outcome of all these options is a list of x,y points with varying spacing and number 
    if option == 'coord':
        geod = pyproj.Geod(ellps='WGS84')
        dlat = 0.3 # at 30km resolution, distance between cells in deg
        dlon = 0.98
        
        if transect_name != '':
            #TODO edit using region_name,region_xybounds, region_coordbounds
            lat = [region_coordbounds[region_name.index(transect_name),1,0],
                   region_coordbounds[region_name.index(transect_name),1,1]]
            lon = [region_coordbounds[region_name.index(transect_name),0,0],
                   region_coordbounds[region_name.index(transect_name),0,1]]
            zmax = region_zbounds[region_name.index(transect_name),1]
            zmin = region_zbounds[region_name.index(transect_name),0]
        else:
            transect_name = (str(int(abs(lat[0]))) + 'S' + 
                             str(int(abs(lon[0]-360))) + 'W-' + 
                             str(int(abs(lat[1]))) + 'S' + 
                             str(int(abs(lon[1]-360))) + 'W' )
            zmax = 0.
            zmin = -9999.
        
        p1,p2 = (lon[0]-360,lat[0]), (lon[1]-360,lat[1])
        frac_along_route = 0.05
        lat_interp = np.zeros((int(1/frac_along_route)))
        lon_interp = np.zeros((int(1/frac_along_route)))
        transect_idx = np.zeros((int(1/frac_along_route)),dtype=int)
            
        for i,frac in enumerate(np.arange(0,1,frac_along_route)):
            lon_interp[i],lat_interp[i] = great_circle.intermediate_point(p1,p2,frac)
            if lon_interp[i] < 0:
                lon_interp[i] += 360
            logical_trans = ( (latcell > (lat_interp[i]-dlat)*deg2rad )    & 
                              (latcell < (lat_interp[i]+dlat)*deg2rad )    &
                              (loncell > (lon_interp[i]-dlon)*deg2rad) & 
                              (loncell < (lon_interp[i]+dlon)*deg2rad) 
                              ) 
            candidate_idx = np.asarray(logical_trans.nonzero(),dtype=int)[0,:]
            distance_from_point = np.zeros((np.shape(candidate_idx)))
            _,_,distance_from_point = geod.inv(
                                      [lon_interp[i] for j in candidate_idx],
                                      [lat_interp[i] for j in candidate_idx],
                                      loncell[candidate_idx]/deg2rad,
                                      latcell[candidate_idx]/deg2rad)
            transect_idx[i] = int(candidate_idx[np.argmin(distance_from_point)])
        _,temp_idx = np.unique(transect_idx,return_index=True)
        cellidx = transect_idx[np.sort(temp_idx)]
        zbool = (zbottom[cellidx] < zmax) & (zbottom[cellidx] > zmin) 
        cellidx = cellidx[zbool]
        edgeidx = []
        dist = np.sqrt( np.square(fmesh.variables['yCell'][cellidx] - 
                                  fmesh.variables['yCell'][cellidx[0]]) + 
                        np.square(fmesh.variables['xCell'][cellidx] - 
                                  fmesh.variables['xCell'][cellidx[0]]) )
        #elif select_cell == 'connecting':
        #    # NOT FUNCTIONAL
        
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

    x_transect = fmesh.variables['xCell'][cellidx]
    y_transect = fmesh.variables['yCell'][cellidx]
    edgesOnCell = np.subtract(fmesh.variables['edgesOnCell'][cellidx,:],1)
    verticesOnCell = np.subtract(fmesh.variables['verticesOnCell'][cellidx,:],1)
    verticesOnEdge = np.zeros((len(cellidx),7,2))
    for i in range(len(cellidx)):
        for j in range(7):
            verticesOnEdge[i,j,:] = (
              fmesh.variables['verticesOnEdge'][edgesOnCell[i,j],:])
    
    # Select edges based on their orientation with respect to the transect line
    edgeidx = [] 
    x0 = xcell[cellidx[0]]
    y0 = ycell[cellidx[0]]
    x1 = xcell[cellidx[-1]]
    y1 = ycell[cellidx[-1]]
    m = (y1 - y0)/(x1 - x0)
    b = (x1*y0 - x0*y1)/(x1 - x0)
    angle = atan(1/m)
    if vartype == 'velocity':
        #transect_angle = angle/deg2rad 
        #if transect_angle>=180:
        #   transect_angle += -360
        #elif transect_angle<-180:
        #   transect_angle += 360
        #if transect_angle < 0:
        #   transect_angle += 180
        #print('transect_angle = ',transect_angle)
        #print('x0,y0 = ',x0,y0)
        #print('x1,y1 = ',x1,y1)
        idx = [] 
        dxy = 1e3
        #dxy = 5e3
        for i,celli in enumerate(cellidx):
            for j,edge in enumerate(edgesOnCell[i,:]):
                angleEdge = fmesh.variables['angleEdge'][edge]#edgesOnCell[i,j]]
                ye = fmesh.variables['yEdge'][edge]
                xe = fmesh.variables['xEdge'][edge]
                if abs(x0-x1) > (y0-y1):#only restrictions on y
                    ym = m*xe + b
                    xlim = xe#xcell[celli]
                    ylim = ym - dxy
                    #print('xe,ym = ',xe,ym)
                else:#only restrictions on x
                    xm = (ye - b)/m
                    xlim = xm - dxy
                    ylim = ye#ycell[celli]
                    #print('xm,ye = ',xm,ye)
                #if i == 1:
                #    print('xcell,ycell = ',xcell[celli],ycell[celli])
                #    print('xe,ye = ',xe,ye)
                #    print('xlim,ylim = ',xlim,ylim)
                if ( (xe <= xlim) and (ye <= ylim) 
                     #and abs(xe-xcell[celli])>1e3 and abs(ye-ycell[celli])>1e3
                   ):
                    #edge not in edgesOnCell[i+1,:]) ):
                    edgeidx.append(edgesOnCell[i,j])
                    idx.append(celli)
                    #edge_angle = angleEdge/deg2rad
                    #if edge_angle>=180:
                    #   edge_angle += -360
                    #elif edge_angle<-180:
                    #   edge_angle += 360
                    #if edge_angle < 0:
                    #   edge_angle += 180
                    #dangle = edge_angle-transect_angle
                    #print('edge_angle = ',edge_angle)
                    #print('dangle = ',dangle)
        cellidx = idx
    if vartype == 'velocity':
        dist = np.sqrt( np.square(fmesh.variables['yEdge'][edgeidx] - 
                                  fmesh.variables['yEdge'][edgeidx[0]]) + 
                        np.square(fmesh.variables['xEdge'][edgeidx] - 
                                  fmesh.variables['xEdge'][edgeidx[0]]) )
    else:
        dist = np.sqrt( np.square(fmesh.variables['yCell'][cellidx] - 
                                  fmesh.variables['yCell'][cellidx[0]]) + 
                        np.square(fmesh.variables['xCell'][cellidx] - 
                                  fmesh.variables['xCell'][cellidx[0]]) )
        
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
    #dist = np.sqrt( np.square(fmesh.variables['yEdge'][edgeidx] - fmesh.variables['yEdge'][edgeidx[0]]) + 
    #                np.square(fmesh.variables['xEdge'][edgeidx] - fmesh.variables['xEdge'][edgeidx[0]]) )
    
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
        idx_scope = pick_from_region(region=scope_name, run=run, plot_map=False)
        zmax_scope = np.multiply(-1.,fmesh.variables['bottomDepth'][idx_scope])
        icemask_scope  = fmesh.variables['landIceMask'][0,idx_scope]
        #zice_scope     = fmesh.variables['landIceDraft'][idx_scope]
        ycell_scope = ycell[idx_scope]
        xcell_scope = xcell[idx_scope]
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cntr1 = plt.scatter(ycell_scope,xcell_scope,
                            s=loc_ptsize[region_name.index(scope_name)], c=zmax_scope)
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
        plt.ylim([region_xybounds[region_name.index(scope_name)][0,0],
                  region_xybounds[region_name.index(scope_name)][0,1]])
        plt.xlim([region_xybounds[region_name.index(scope_name)][1,0],
                  region_xybounds[region_name.index(scope_name)][1,1]])
        cbar = plt.colorbar(cntr1)
        #cbar = plt.colorbar(sc)
        cbar.set_label(r'Depth (m)')    
        print('save ','bathy_' + transect_name)
        plt.savefig(savepath + 'bathy_' + transect_name)
        plt.close()
    
    return cellidx, edgeidx, dist, angle

def pick_from_region(region='frisEAcoast',run = 'ISMF',
                     land_ice_mask = False, plot_map=False,
                     overwrite=False, savepath=savepath_nersc):
    
    fmesh = netCDF4.Dataset(meshpath[runname.index(run)])

    # import variables from file
    latCell = fmesh.variables['latCell'][:]
    lonCell = fmesh.variables['lonCell'][:]
    xCell   = fmesh.variables['xCell'][:]
    yCell   = fmesh.variables['yCell'][:]
    zmax    = np.multiply(-1,fmesh.variables['bottomDepth'][:])
    idx_bool = np.ones_like(xCell, dtype=bool)
    if region_xybounds[region_name.index(region)] is not None:
        idx_bool = (idx_bool & 
                    (xCell < region_xybounds[region_name.index(region),0,1]) & 
                    (xCell > region_xybounds[region_name.index(region),0,0]) & 
                    (yCell < region_xybounds[region_name.index(region),1,1]) & 
                    (yCell > region_xybounds[region_name.index(region),1,0]))
    if region_coordbounds[region_name.index(region)] is not None:
        idx_bool = (idx_bool & 
                    (latCell < region_coordbounds[region_name.index(region),1,1]*deg2rad) & 
                    (latCell > region_coordbounds[region_name.index(region),1,0]*deg2rad) & 
                    (lonCell < region_coordbounds[region_name.index(region),0,1]*deg2rad) & 
                    (lonCell > region_coordbounds[region_name.index(region),0,0]*deg2rad))
    if region_zbounds[region_name.index(region)] is not None:
        idx_bool = (idx_bool & 
                    (zmax < region_zbounds[region_name.index(region),1]) & 
                    (zmax > region_zbounds[region_name.index(region),0]))
    if land_ice_mask:
        landIceMask = fmesh.variables['landIceMask'][:]
        idx_bool = (idx_bool & (landIceMask == 1))
    cellidx = np.asarray(idx_bool.nonzero(),dtype=int)[0,:]
    return cellidx
