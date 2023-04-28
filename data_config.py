#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:56:52 2019

@author: cbegeman
"""

runname = ['CGM-UIB', 'CGM-DIB', 'VGM-DIB']
savepath = '/lcrc/group/e3sm/ac.cbegeman/scratch/e3sm_fris_analysis'
obspath = '/lcrc/group/acme/diagnostics/observations'
runpath = ['/lcrc/group/acme/ac.mhoffman/acme_scratch/anvil/20191003.GMPAS-IAF-ISMF.T62_oEC60to30v3wLI.cori-knl/archive/ocn/hist',
           '/lcrc/group/acme/ac.mhoffman/acme_scratch/anvil/20210730.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.DIBbugfix.anvil/run',
           '/lcrc/group/acme/ac.mhoffman/acme_scratch/anvil/20210901.GMPAS-DIB-IAF-ISMF.T62_oEC60to30v3wLI.VGM.DIBbugfix.anvil/run']
meshpath = ['/lcrc/group/acme/data/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc',
            '/lcrc/group/acme/data/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc',
            '/lcrc/group/acme/data/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc']
