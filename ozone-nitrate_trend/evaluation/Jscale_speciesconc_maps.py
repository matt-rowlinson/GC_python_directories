#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=map_plots
#SBATCH --ntasks=1
#SBATCH --mem=45gb
#SBATCH --partition=nodes
#SBATCH --time=00:25:00
#SBATCH --output=Logs/map_plots_%A.log
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


def get_model_as_dataframe(rundirs, variable, version="13.1.2", collection="SpeciesConc", year="2017", x=-24.9, y=16.9, z=0):  
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

site_dict = {
        "CVAO" : {
            "longitude" : -24.9,
            "latitude"  : 16.9,
            "altitude"  : 0 },
        "Mace_Head" : {
            "longitude" : -24.9,
            "latitude"  : 16.9,
            "altitude"  : 0 }           
        }

###----------------------------------------MAIN SCRIPT----------------------------------------------------------------------##

def main():
    #variables=#O3, NIT, NITs, SO2, SO4,'ethane','propane','ALK4','benzene','toluene','EOH'] 
    variables=['CO']
    for variable in variables:
        print(variable)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir='nitrate_photol_Jscale-25', variable=variable, 
                                                  version='13.1.2', 
                                                  year='')    
        
        control_df = var[:1461]
        dev_df     = var[-1461:]

        control_df = np.mean(control_df[:,0,:,:],0)
        dev_df     = np.mean(dev_df[:,0,:,:],0)

        fig,(ax1, ax2, ax3)= plt.subplots(1,3,figsize=(12,3))
        X,Y=np.meshgrid(lon,lat)
        m1=rp.get_basemap(ax=ax1, lines=False, freq=60.)
        m2=rp.get_basemap(ax=ax2, lines=False, freq=60.)
        m3=rp.get_basemap(ax=ax3, lines=False, freq=60.)

        MAX=np.nanmax(control_df)

        diff = - (100 - ( dev_df / control_df ) * 100 )
        Min, Max = rp.get_abs_max_min(diff)
        if Max > 1000:
            Min=-1000 ; Max=1000
        
        im=ax1.pcolormesh(X,Y, control_df,  cmap=cmap, vmax=MAX)
        im2=ax2.pcolormesh(X,Y, dev_df,  cmap=cmap, vmax=MAX)
        im3=ax3.pcolormesh(X,Y, diff,  cmap='bwr', vmin=Min, vmax=Max)#, vmin=Min)

        ax1.set_title('a) 1980-1983')
        ax2.set_title('b) 2017-2020')
        ax3.set_title('b / a')

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.2, top=.9)
        cbar_ax = fig.add_axes([0.065, 0.15, 0.552, 0.07])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel(f'{variable} ()')

        cbar_ax = fig.add_axes([0.69,0.15,0.285, 0.07])
        fig.colorbar(im3, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel(f'%')#delta {variable}')

        plt.suptitle(f'SpeciesConc_{variable}')
        plt.savefig(f'plots/SpeciesConc_map_1980_2020_{variable}.png')
        plt.close()

if __name__ == "__main__":
    main()
