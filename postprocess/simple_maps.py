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
import glob
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

variable='CO'
out='COconc'

filen='/mnt/lustre/users/mjr583/GC/GCHP/rundirs/test_run/sg_CV_3hour/OutputDir/'
for infile in glob.glob(filen+'*SpeciesConc*0030*nc*'):
    print(infile)

from netCDF4 import Dataset
fh=Dataset(infile)
O3=np.squeeze(fh.variables['SpeciesConc_O3'][:])
x=fh.variables['Xdim'][:]
y=fh.variables['Ydim'][:]
print(x)
print(y)

O3=O3[0]
X,Y = np.meshgrid(x,y)

fig,((ax1, ax2, ax3),(ax4,ax5,ax6))= plt.subplots(2,3,figsize=(12,7))
axes=[ax1, ax2, ax3,ax4,ax5,ax6]

for n,ax in enumerate(axes):
    #m=rp.get_basemap(lllon=-70, urlon=0., lllat=-15., urlat=30.,ax=ax1)
    ax.pcolormesh(X,Y,O3[n],cmap=cmap)

plt.savefig('./XXX.png')
sys.exit()
