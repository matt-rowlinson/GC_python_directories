#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=map_plots
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --partition=interactive
#SBATCH --time=00:15:00
#SBATCH --output=Logs/9_panel_map_plots_%A.log
import os
import sys
import numpy as np
import xarray as xr
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

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', years=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    if years:
        f_list=[]
        for y in years:
            f_list.append( find_file_list(path, [f'SpeciesConc.{y}']) )
            ds = xr.open_mfdataset( f_list, combine='by_coords' )
    else:
        ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    return ds


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


###----------------------------------------MAIN SCRIPT----------------------------------------------------------------------##

def main():
    #variables=#O3, NIT, NITs, SO2, SO4,'ethane','propane','ALK4','benzene','toluene','EOH'] 
    variables=['O3']
    Abs=False ; suff=''
    rundir='trimmed_nitrate_photol_Jscale-25'
    
    fig,((ax1, ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8),(ax9,ax10))= plt.subplots(5,2,figsize=(9,12))
    for variable in variables:
        print(variable)
        lon = get_data_as_xr(rundir, year='198001')['lon']#.mean(dim='time')[0]
        lat = get_data_as_xr(rundir, year='198001')['lat']#.mean(dim='time')[0]

        control = get_data_as_xr(rundir, year='1980')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev00 = get_data_as_xr(rundir, year='2000')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev05 = get_data_as_xr(rundir, year='2005')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev10 = get_data_as_xr(rundir, year='2010')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev17 = get_data_as_xr(rundir, year='2017')[f'SpeciesConc_{variable}'].mean(dim='time')[0]

        print( control.shape )
        print( dev10.shape )


         

        X,Y=np.meshgrid(lon,lat)
        m1=rp.get_basemap(ax=ax1,  lines=False, freq=60.)
        m2=rp.get_basemap(ax=ax3,  lines=False, freq=60.)
        m3=rp.get_basemap(ax=ax4,  lines=False, freq=60.)
        m4=rp.get_basemap(ax=ax5,  lines=False, freq=60.)
        m5=rp.get_basemap(ax=ax6,  lines=False, freq=60.)
        m6=rp.get_basemap(ax=ax7,  lines=False, freq=60.)
        m7=rp.get_basemap(ax=ax8,  lines=False, freq=60.)
        m8=rp.get_basemap(ax=ax9,  lines=False, freq=60.)
        m9=rp.get_basemap(ax=ax10, lines=False, freq=60.)

        MAX=np.nanmax([control.max(),dev00.max(),dev05.max(),dev10.max(),dev17.max() ] )
        
        ax1.set_title('1980')
        ax3.set_title('2000')
        ax5.set_title('2005')
        ax7.set_title('2010') 
        ax9.set_title('2017') 

        if Abs:
            diff90 = dev00 - control
            diff00 = dev05 - control
            diff10 = dev10 - control
            diff20 = dev17 - control
            ax4.set_title( '2000 / 1980')
            ax6.set_title( '2005 / 1980')
            ax8.set_title( '2010 / 1980')
            ax10.set_title('2017 / 1980')

            suff='_abs'
            clabel=f"delta {variable} ({d[variable]['unit']})"

        else:
            diff90 = - (100 - ( dev00 / control ) * 100 )
            diff00 = - (100 - ( dev05 / control ) * 100 )
            diff10 = - (100 - ( dev10 / control ) * 100 )
            diff20 = - (100 - ( dev17 / control ) * 100 )
            ax4.set_title( '2000 / 1980')
            ax6.set_title( '2005 / 1980')
            ax8.set_title( '2010 / 1980')
            ax10.set_title('2017 / 1980')
            clabel=f"delta {variable} (%)"

        Min90, Max90 = rp.get_abs_max_min(diff90)
        Min00, Max00 = rp.get_abs_max_min(diff00)
        Min10, Max10 = rp.get_abs_max_min(diff10)
        Min20, Max20 = rp.get_abs_max_min(diff20)
        Max = np.max([Max90, Max00, Max10, Max20])
        Min = -Max
        if Max > 1000:
            Min=-1000 ; Max=1000
        
        im=ax1.pcolormesh(X,Y, control,  cmap=cmap, vmax=MAX)
        im=ax3.pcolormesh(X,Y, dev00,  cmap=cmap, vmax=MAX)
        im=ax5.pcolormesh(X,Y, dev05,  cmap=cmap, vmax=MAX)
        im=ax7.pcolormesh(X,Y, dev10,  cmap=cmap, vmax=MAX)
        im=ax9.pcolormesh(X,Y, dev17,  cmap=cmap, vmax=MAX)

        im3=ax4.pcolormesh(X,Y, diff90,  cmap='bwr', vmin=Min, vmax=Max)
        im3=ax6.pcolormesh(X,Y, diff00,  cmap='bwr', vmin=Min, vmax=Max)
        im3=ax8.pcolormesh(X,Y, diff10,  cmap='bwr', vmin=Min, vmax=Max)
        im3=ax10.pcolormesh(X,Y, diff20,  cmap='bwr', vmin=Min, vmax=Max)

        plt.delaxes(ax2)


        plt.tight_layout() 
        plt.subplots_adjust(bottom=0.1, top=0.94)
        cbar_ax = fig.add_axes([0.05, 0.07, 0.4, 0.014])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel(f'{variable} ({d[variable]["unit"]})')

        cbar_ax = fig.add_axes([0.55,0.07,0.4, 0.014])
        fig.colorbar(im3, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel(clabel)

        plt.suptitle(f'SpeciesConc_{variable}', weight='bold')
        plt.savefig(f'plots/{rundir}_decadal_maps_{variable}{suff}.png', dpi=400)
        plt.close()

if __name__ == "__main__":
    main()
