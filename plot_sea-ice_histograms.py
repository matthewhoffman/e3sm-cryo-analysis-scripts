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

run = ['ISMF','ISMF-noEAIS','ISMF-noDIB','ISMF-3dGM']
runlabel = ['CTRL','noEAmelt','uniformIB','modGM']
savepath = '/Users/cbegeman/Data/weddell_analysis/'
filename1 = 'ISMF_ISMF-noEAIS_ISMF-3dGM_ISMF-noDIB_sea_ice_fw_flux_frisEAcoast_t010-180.txt'
ln_color = ['k','#1b9e77','#7570b3','#d95f02']
plot_filename = ''
s_yr = 3600.*24.*365.25
fig = plt.figure()
ax = fig.add_subplot(111)
df = pandas.read_csv(savepath + filename1)
_t = df['decyear'][:]
for k,run_incr in enumerate(run):
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
fig.savefig('/Users/cbegeman/Data/weddell_analysis/G-ISMF_sea_ice_fw_flux_pdf.pdf')
