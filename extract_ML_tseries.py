#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
'''

import sys
import numpy as np
import weddell_mod as wed
from plot_config import savepath_anvil

wed.tseries1(['ISMF-3dGM'], ['tMLD','bMLD','BLD','surfBuoy'],
             year_range=[30,200], 
             lat=61.5, lon=360-58,
             print_to_file=True, create_figure=False,
             savepath=savepath_anvil, overwrite=True)

