#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=15map_plots
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --partition=interactive
#SBATCH --time=00:15:00
#SBATCH --output=Logs/15_panel_map_plots_%A.log
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
    Abs=True ; suff='' ; scale='100'
    
    fig,((ax1, ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(ax10,ax11,ax12),(ax13,ax14,ax15))= plt.subplots(5,3,figsize=(12,12))
    axes=[ax1, ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15]
    titles=['Control 1980',f'J-{scale} 1980', f'J-{scale} - control',
            'Control 2000',f'J-{scale} 2000', f'J-{scale} - control',
            'Control 2005',f'J-{scale} 2005', f'J-{scale} - control',
            'Control 2010',f'J-{scale} 2010', f'J-{scale} - control',
            'Control 2017',f'J-{scale} 2017', f'J-{scale} - control' ]

    for variable in variables:
        print(variable)

        rundir='trimmed_nitrate_photol_control'
        lon = get_data_as_xr(rundir, year='198001')['lon']#.mean(dim='time')[0]
        lat = get_data_as_xr(rundir, year='198001')['lat']#.mean(dim='time')[0]
        dev80 = get_data_as_xr(rundir, year='1980')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev00 = get_data_as_xr(  rundir, year='2000')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev05 = get_data_as_xr(  rundir, year='2005')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev10 = get_data_as_xr(  rundir, year='2010')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        dev17 = get_data_as_xr(  rundir, year='2017')[f'SpeciesConc_{variable}'].mean(dim='time')[0]

        
        rundir=f'trimmed_nitrate_photol_Jscale-{scale}'
        Jdev80 = get_data_as_xr(rundir, year='1980')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        Jdev00 = get_data_as_xr(rundir, year='2000')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        Jdev05 = get_data_as_xr(rundir, year='2005')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        Jdev10 = get_data_as_xr(rundir, year='2010')[f'SpeciesConc_{variable}'].mean(dim='time')[0]
        Jdev17 = get_data_as_xr(rundir, year='2017')[f'SpeciesConc_{variable}'].mean(dim='time')[0]


        X,Y=np.meshgrid(lon,lat)
        ax1.set_title('Base')
        ax2.set_title(f'Jscale-{scale}')

        for n, ax in enumerate(axes):
            m=rp.get_basemap(ax=ax,  lines=False, freq=60.)
            #ax.set_title(titles[n])

        MAX=np.nanmax([dev80.max(),dev00.max(),dev05.max(),dev10.max(),dev17.max(),
                       Jdev80.max(),Jdev00.max(),Jdev05.max(),Jdev10.max(),Jdev17.max() ] )
        

        if Abs:
            diff80 = Jdev80 - dev80
            diff00 = Jdev00 - dev00
            diff05 = Jdev05 - dev05
            diff10 = Jdev10 - dev10
            diff17 = Jdev17 - dev17
            suff='_abs'
            clabel=f"delta {variable} ({d[variable]['unit']})"
            ax3.set_title( 'Jscale-{scale} - Base')

        else:
            diff80 = - (100 - ( Jdev80 / dev80 ) * 100 )
            diff00 = - (100 - ( Jdev00 / dev00 ) * 100 )
            diff05 = - (100 - ( Jdev05 / dev05 ) * 100 )
            diff10 = - (100 - ( Jdev10 / dev10 ) * 100 )
            diff17 = - (100 - ( Jdev17 / dev17 ) * 100 )
            ax3.set_title( 'Jscale-{scale} / Base')
            clabel=f"delta {variable} (%)"

        Min80, Max80 = rp.get_abs_max_min(diff80)
        Min00, Max00 = rp.get_abs_max_min(diff00)
        Min05, Max05 = rp.get_abs_max_min(diff05)
        Min10, Max10 = rp.get_abs_max_min(diff10)
        Min20, Max20 = rp.get_abs_max_min(diff17)
        Max = np.max([Max80, Max00, Max05, Max10, Max20])
        Min = -Max
        if Max > 1000:
            Min=-1000 ; Max=1000
        
        im=ax1.pcolormesh(X,Y, dev80,  cmap=cmap, vmax=MAX)
        im=ax4.pcolormesh(X,Y, dev00,  cmap=cmap, vmax=MAX)
        im=ax7.pcolormesh(X,Y, dev05,  cmap=cmap, vmax=MAX)
        im=ax10.pcolormesh(X,Y, dev10,  cmap=cmap, vmax=MAX)
        im=ax13.pcolormesh(X,Y, dev17,  cmap=cmap, vmax=MAX)
        
        im=ax2.pcolormesh(X,Y, Jdev80,  cmap=cmap, vmax=MAX)
        im=ax5.pcolormesh(X,Y, Jdev00,  cmap=cmap, vmax=MAX)
        im=ax8.pcolormesh(X,Y, Jdev05,  cmap=cmap, vmax=MAX)
        im=ax11.pcolormesh(X,Y, Jdev10,  cmap=cmap, vmax=MAX)
        im=ax14.pcolormesh(X,Y, Jdev17,  cmap=cmap, vmax=MAX)

        im3=ax3.pcolormesh(X,Y,  diff80,  cmap='bwr', vmin=Min, vmax=Max)
        im3=ax6.pcolormesh(X,Y,  diff00,  cmap='bwr', vmin=Min, vmax=Max)
        im3=ax9.pcolormesh(X,Y,  diff05,  cmap='bwr', vmin=Min, vmax=Max)
        im3=ax12.pcolormesh(X,Y, diff10,  cmap='bwr', vmin=Min, vmax=Max)
        im3=ax15.pcolormesh(X,Y, diff17,  cmap='bwr', vmin=Min, vmax=Max)

        plt.tight_layout()
        plt.subplots_adjust(left=.12, hspace=.15)
        #ax.text(-1., 0., '1980',  rotation=90, fontsize=200)
        #ax.text(.1, 0.5, '2000',  rotation=90, fontsize=200)
        #ax.text(.1, 0.5, '2005',  rotation=90, fontsize=200)
        #ax.text(.1, 0.5, '2010',  rotation=90, fontsize=200)
        #ax.text(.1, 0.5, '2017',  rotation=90, fontsize=200)

        plt.subplots_adjust(bottom=0.1, top=0.94)
        cbar_ax = fig.add_axes([0.05, 0.07, 0.4, 0.014])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel(f'{variable} ({d[variable]["unit"]})')

        cbar_ax = fig.add_axes([0.55,0.07,0.4, 0.014])
        fig.colorbar(im3, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel(clabel)

        plt.suptitle(f'SpeciesConc_{variable}', weight='bold')
        plt.savefig(f'plots/{rundir}_15decadal_maps_{variable}{suff}.png', dpi=400)
        plt.close()

if __name__ == "__main__":
    main()
