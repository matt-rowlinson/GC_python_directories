#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=pcolor_plot
#SBATCH --ntasks=1
#SBATCH --mem=15gb
#SBATCH --partition=nodes
#SBATCH --time=00:11:30
#SBATCH --output=LOGS/pcolormesh_%a.log
#SBATCH --array=1-2
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import matplotlib.colors as mcolors
from CVAO_dict import CVAO_dict as d

current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd


from GC_tools import Bry as D
variable='Bry'
a, lat, lon, lev, atime = GC.get_var_group(D, rundir='tropchem_merra_4x5', version='12.9.3', year='2016')
b, lat,lon,lev,btime = GC.get_var_group(D, rundir='fullchem_4x5_LVOCfalse', version='GEOS-Chem',year='2016')

a=np.mean(a[:,0,:,:],0)
b=np.mean(b[:,0,:,:],0)

rp.difference_map(a,b, variable, lon, lat, labels=['v12.9.3','v13.0.0','v12 - v13'])
rp.difference_map(a,b, variable, lon, lat, labels=['v12.9.3','v13.0.0','v12 - v13'], pc=True)
