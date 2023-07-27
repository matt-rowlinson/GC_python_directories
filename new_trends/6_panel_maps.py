#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=map_plots
#SBATCH --ntasks=1
#SBATCH --mem=45gb
#SBATCH --partition=nodes
#SBATCH --time=00:45:00
#SBATCH --output=Logs/9_panel_map_plots_%A.log
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import matplotlib.colors as mcolors
from CVAO_dict import CVAO_dict as d
plt.style.use('seaborn-darkgrid')

current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd


def get_model_as_dataframe(rundirs, variable, version="13.1.2", collection="SpeciesConc", year="2018", x=-24.9, y=16.9, z=0):  
    dfs=[]
    for a in rundirs:
        print(a)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, collection=collection, version='13.1.2', year=year)
        #var = var * 1e12
        a_time=[] 
        idy = rp.find_nearest(lat, y)
        idx = rp.find_nearest(lon, x)
        for t in range(len(atime)):
            v=var[t,z,idy,idx]
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    cv_df=pd.concat(dfs, axis=1, ignore_index=False)
    return dfs


###----------------------------------------MAIN SCRIPT----------------------------------------------------------------------##

def main():
    variables=['NIT','NITs','NITD1','NITD2','NITD3','NITD4']
    
    fig,((ax1, ax2),(ax3,ax4),(ax5,ax6))= plt.subplots(3,2,figsize=(16,12))
    axes=[ax1, ax2,ax3,ax4,ax5,ax6]
    maps=[]
    for n, variable in enumerate(variables):
        print(variable)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir='nitrate_photol_control', variable=variable, 
                                                  version='13.1.2_DUST', 
                                                  year='2018')    
        var = np.mean(var[:,0,:,:],0) * 1e12
        maps.append(var)

    for n, ax in enumerate(axes):
        X,Y=np.meshgrid(lon,lat)
        m1=rp.get_basemap(ax=ax,  lines=False, freq=60.)

        Max = np.max(maps)
        mean = str(np.round( np.max( maps[n] ), 2))
        cv_mean = str(np.round( np.mean( maps[n][31,27] ), 2))
        ax.set_title(f'{variables[n]} (mean={mean} pptV)', fontsize=18)
        im=ax.pcolormesh(X,Y, maps[n],  cmap=cmap, norm=LogNorm(),vmin=0.1, vmax=Max)
        ax.scatter(lon[31],lat[27], marker='x', color='k')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='horizontal')

    plt.tight_layout()
    plt.savefig(f'plots/NIT_global_surface.png', dpi=400)
    plt.close()

if __name__ == "__main__":
    main()
