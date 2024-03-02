#!/usr/bin/env python
import numpy as np

import weddell_mod as wed
from data_config import runname
import mpas_analysis.shared.plot.inset as inset

#runs=['CGM-DIB', 'VGM-DIB']; years = [191]; year_range=np.arange(191,201,1)
#runs=['CGM-UIB',]; years = [91]; year_range=np.arange(91,101,1)

#runs=['CGM-DIB', 'VGM-DIB']; years = [191]; year_range=np.arange(191,192,1)
runs=['CGM-UIB',]; years = [91]; year_range=np.arange(91,92,1)

months = [1]
rho_levels = np.arange(1027.4,1028.0,0.1)
sigma_levels = np.arange(31.3,32.7,0.05)
sigma_levels_label = np.linspace(32.0, 32.7, 8)
print(sigma_levels_label)

for run in runs:
    fig = wed.transect('coord', years, months, ['T'], 
                 run=run, transect_name='trough_crossshelf',
                 plot_transect_on_map=True,
                 year_range=year_range,
                 month_range=np.arange(1,13,1),
                 ops=['time_mean'], overwrite=True, 
                 var_contour='sigma1', cntr_levels=sigma_levels,
                 cntr_label_levels=sigma_levels_label,
                 zlim=[5,-1e3], figure_format='pdf')

    