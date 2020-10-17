#!/usr/bin/env python

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm
from math import pi
import weddell_mod as wed
#from transects import transect

wed.plot_surf_var('rho',98,1,level=1027.6,run='ISMF')
