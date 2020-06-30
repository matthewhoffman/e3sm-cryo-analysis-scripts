#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

import sys
import os
import netCDF4
import numpy as np

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
    return zmid,zbottom,ztop
