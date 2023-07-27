#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=diff_diff
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --partition=interactive
#SBATCH --time=00:15:00
#SBATCH --output=Logs/diff_diff_%A.log
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
    Abs=True ; suff='' ; scale='25'
    
    fig,((ax1,ax2,ax3,ax4))= plt.subplots(1,4,figsize=(11,3))
    axes=[ax1, ax2,ax3,ax4]
    titles=['2000 - 1980','2005 - 1980','2010 - 1980','2017 - 1980']

    for variable in variables:
        print(variable)

        rundir='trimmed_nitrate_photol_control'
        lon = get_data_as_xr(rundir, year='198001')['lon']#.mean(dim='time')[0]
        lat = get_data_as_xr(rundir, year='198001')['lat']#.mean(dim='time')[0]
        dev80 = get_data_as_xr(rundir, year='1980')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        dev00 = get_data_as_xr(  rundir, year='2000')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        dev05 = get_data_as_xr(  rundir, year='2005')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        dev10 = get_data_as_xr(  rundir, year='2010')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        dev17 = get_data_as_xr(  rundir, year='2017')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9

        rundir=f'trimmed_nitrate_photol_Jscale-{scale}'
        Jdev80 = get_data_as_xr(rundir, year='1980')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        Jdev00 = get_data_as_xr(rundir, year='2000')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        Jdev05 = get_data_as_xr(rundir, year='2005')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        Jdev10 = get_data_as_xr(rundir, year='2010')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9
        Jdev17 = get_data_as_xr(rundir, year='2017')[f'SpeciesConc_{variable}'].mean(dim='time')[0] * 1e9

        X,Y=np.meshgrid(lon,lat)
        for n, ax in enumerate(axes):
            m=rp.get_basemap(ax=ax,  lines=False, freq=60.)
            ax.set_title(titles[n])

        if Abs:
            diff80 = Jdev80 - dev80
            diff00 = Jdev00 - dev00
            diff05 = Jdev05 - dev05
            diff10 = Jdev10 - dev10
            diff17 = Jdev17 - dev17
            suff='_abs'
            clabel=f"delta {variable} ({d[variable]['unit']})"

        else:
            diff80 = - (100 - ( Jdev80 / dev80 ) * 100 )
            diff00 = - (100 - ( Jdev00 / dev00 ) * 100 )
            diff05 = - (100 - ( Jdev05 / dev05 ) * 100 )
            diff10 = - (100 - ( Jdev10 / dev10 ) * 100 )
            diff17 = - (100 - ( Jdev17 / dev17 ) * 100 )
            clabel=f"delta {variable} (%)"

        to_plot = [ diff00 - diff80, diff05 - diff80, diff10 - diff80, diff17 - diff80,]
        
        Min00, Max00 = rp.get_abs_max_min(diff00 - diff80)
        Min05, Max05 = rp.get_abs_max_min(diff05 - diff80)
        Min10, Max10 = rp.get_abs_max_min(diff10 - diff80)
        Min20, Max20 = rp.get_abs_max_min(diff17 - diff80)
        Max = np.max([Max00, Max05, Max10, Max20])
        Min = -Max
        if Max > 1000:
            Min=-100 ; Max=100

        for n, ax in enumerate(axes):
            im=ax.pcolormesh(X,Y,  to_plot[n],  cmap='bwr', vmin=Min, vmax=Max)
            print( 'Mean :', to_plot[n].values.mean() )

        plt.tight_layout()
        #plt.subplots_adjust(bottom=0.1 )
        cbar_ax = fig.add_axes([0.15, 0.2, 0.7, 0.06])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar_ax.set_xlabel(clabel)

        plt.suptitle(f'$O_3$ sensitivity to Jscale-{scale}', weight='bold')
        #plt.tight_layout()
        plt.savefig(f'plots/diff_diff_{rundir}_{variable}{suff}.png', dpi=400)
        plt.close()

if __name__ == "__main__":
    main()
