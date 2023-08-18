#!/usr/bin/env python
import numpy as np

import weddell_mod as wed
from data_config import runname

years = [191]
months = [1]
rho_levels = np.arange(1027.4,1028.0,0.1)
sigma_levels = np.arange(32.0,32.7,0.05)

for run in ['CGM-DIB', 'VGM-DIB']:
    wed.transect('coord', years, months, ['S'], 
                 run=run, transect_name='trough_crossshelf',
                 year_range=np.arange(191,201,1),
                 month_range=np.arange(1,13,1),
                 ops=['time_mean'], overwrite=True, 
                 var_contour='sigma1', cntr_levels=sigma_levels,
                 zlim=[5,-1e3])
