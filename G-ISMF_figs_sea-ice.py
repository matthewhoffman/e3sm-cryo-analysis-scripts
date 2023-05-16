#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 15:31:15 2021


@author: cbegeman
"""
import matplotlib
import matplotlib.pyplot as plt
import pandas
import numpy as np
from plot_config import savepath_anvil

run = ['ISMF','ISMF-noDIB','ISMF-3dGM']
runlabel = ['CTRL','uniformIB','modGM']
filename1 = 'ISMF-noDIB_data/ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_sea_ice_fw_flux_frisEAcoast_t010-180.txt'
filename2 = 'ISMF_ISMF-3dGM_sea_ice_fw_flux_frisEAcoast_t001-200.txt'
ln_color = ['k','#7570b3','#d95f02']
plot_filename = ''
s_yr = 3600.*24.*365.25
fig = plt.figure()
ax = fig.add_subplot(111)
filename=[filename2,filename1,filename2]
for k,run_incr in enumerate(run):
    df = pandas.read_csv(savepath_anvil + filename[k])
    t = df['decyear'][:]
    var1 = df[run_incr+'_sea_ice_fw_flux'][:]
    counts, bins = np.histogram(var1,bins=20,
                range=(np.min(var1)*1.2,np.max(var1)*1.1))
    bin_width = bins[1:]-bins[:-1]
    bin_center = 0.5*(bins[1:]+bins[:-1])
    bin_height = (counts/np.sum(counts))/bin_width
    #print(np.sum(counts))
    #print(np.sum(bin_height*bin_width))
    #plt.hist(curl,color=ln_color[k],facealpha=1)
    plt.plot(bin_center,bin_height,ln_color[k],label = runlabel[k])
    del bin_center, bin_height
del df,t
plt.legend()
plt.xlabel(r'Sea ice freshwater flux ($kg \: s^{-1}$)')
plt.ylabel('Probability density')
#plt.ylim([-1.2*ymax,1.2*ymax])
plt.show()
fig.savefig(savepath_anvil+'G-ISMF_sea_ice_fw_flux_pdf.png')
