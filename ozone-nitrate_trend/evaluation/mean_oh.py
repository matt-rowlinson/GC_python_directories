#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/jscale_time_plots_%A.log
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
from CVAO_dict import CVAO_dict as d
import RowPy as rp
import argparse
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
plt.style.use('seaborn-darkgrid')
import constants as const

import matplotlib.colors as mcolors
current_dir = os.path.dirname(__file__)
rgb_WhGrYlRd = np.genfromtxt('/users/mjr583/python_lib/colormaps/WhGrYlRd.txt',
                                     delimiter=' ')
WhGrYlRd = mcolors.ListedColormap(rgb_WhGrYlRd/255.0)
cmap=WhGrYlRd

def get_arguments():
    parser = argparse.ArgumentParser(description="Parse arguments to pass to GC processing scripts")
    parser.add_argument("-r", "--rundir", type=str, 
                        help='Name of desired GC rundir')
    parser.add_argument("-v", "--var", type=list,
                        default=["O3"],
                        help="Name of GC variable")
    parser.add_argument("-V", "--version", type=str,
                        default='13.1.2',
                        help="Version of GEOS-Chem")
    parser.add_argument("-s", "--site", type=str,
                        default="CVO",
                        help="GAW site of interest")
    args=parser.parse_args()
    return args.rundir, args.var, args.version, args.site

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'Metrics.{year}']), combine='by_coords' )
    if variable:
        ds = ds[f'SpeciesConc_{variable}']#.isel(lev=0) * float(d[variable]['scale'])
    return ds

def plot_max_min(plottable):
    ax_Max=[]
    for i in plottable:
        Min, Max = rp.get_abs_max_min(i)
        ax_Max.append( Max.values )
    return ax_Max

def plot(f, axes, plottable, lat, lon, c='', labels=False,
        cbar_label=[],sname="4_panel_plot.png",
        panel_labels=['(a)','(b)','(c)','(d)']):

    X, Y = np.meshgrid( lon, lat )
    ax_Max = plot_max_min(plottable)
    for n in range(len(plottable)):
        m = rp.get_basemap(freq=60, ax=axes[n])
        if n <2:
            im = m.pcolormesh(X, Y, plottable[n], cmap=cmap, vmax=np.max( [ax_Max[0], ax_Max[1] ] ))
        else:
            im = m.pcolormesh(X, Y, plottable[n], cmap='bwr', vmax=ax_Max[n], vmin=-(ax_Max[n]),)
        if labels:
            axes[n].set_title(f'{panel_labels[n]} {labels[n]}', fontsize=14 )
        
        divider = make_axes_locatable(axes[n])
        cax = divider.append_axes('bottom', size='5%', pad=0.25)
        cbar = f.colorbar(im, cax=cax, orientation='horizontal')
        cbar.ax.set_xlabel(f'{cbar_label[n]}',fontsize=12)

    plt.tight_layout()
    plt.savefig( f"plots/{sname}" )
    plt.close()
    return

def oh_interhemispheric_ratio(oh, airmass, lat, mw_oh_kg):
    SHmean_oh = np.nansum(oh.values[0,:23,:]) / np.nansum(airmass.values[0,:23,:])
    SHmean_oh *= (const.AVOGADRO / (mw_oh_kg * 1.0e6)) * 1.0e-5

    NHmean_oh = np.nansum(oh.values[0,23:,:]) / np.nansum(airmass.values[0,23:,:])
    NHmean_oh *= (const.AVOGADRO / (mw_oh_kg * 1.0e6)) * 1.0e-5

    nhsh = NHmean_oh / SHmean_oh
    return nhsh


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    global d, rundirs, v
    rundir, variables, version, site = get_arguments()
    from CVAO_dict import CVAO_dict as d
    sys.path.append( '/users/mjr583/cvao' )
    from lowess_smoother import loess

    rundirs = ['nitrate_photol_2010_control','nitrate_photol_2015_NITS_scale-25']
            #'nitrate_photol_2015_NIT_scale-25','nitrate_photol_2015_NITS_scale-25'     ]
    mw_oh_kg = 17.008 * 1.0e-3
    f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
    axes=[ax1,ax2,ax3,ax4]
    mean_oh_=[] ; diff=[] ; nhsh=[] ;  add_plottable=[]
    for n, rundir in enumerate(rundirs):
        print( n, rundir )
        ds0 = get_data_as_xr(rundir, year='201607')
        lat = ds0['lat']
        lon = ds0['lon']

        oh     = ds0['OHwgtByAirMassColumnFull']
        airmass  = ds0['AirMassColumnFull']

        nhsh.append(np.round( oh_interhemispheric_ratio(oh, airmass, lat, mw_oh_kg),2))

        mean_oh = np.nansum(oh.values) / np.nansum(airmass.values)
        mean_oh *= (const.AVOGADRO / (mw_oh_kg * 1.0e6)) * 1.0e-5
        
        mean_oh_.append(mean_oh)
        oh_2d = oh[0] / airmass[0]
        add_plottable.append( oh[0] )
        if n > 0:
            diff.append(- (100 - ( mean_oh_[1] / mean_oh_[0] * 100 ) ) )
        
    
    add_plottable.append( add_plottable[1] - add_plottable[0] ) 
    add_plottable.append( - ( 100 - ( add_plottable[1] / add_plottable[0] * 100 )))

    plot( f, axes, add_plottable, lat, lon, c=cmap, 
                cbar_label=['Column OH (kg)', 'Column OH (kg)', 'Delta column OH (kg)', '%'],
                labels=[f'Base (NH/SH={nhsh[0]})', f'{rundir.replace("nitrate_photol_2015_","")} (NH/SH={nhsh[1]})', 
                    'Jscale - Base',f'Jscale / Base (mean={np.round(diff[0],2)}%)'],
                sname=f'mean_oh_{rundir.replace("nitrate_photol_2015_","")}')
        

    print( mean_oh_, diff, nhsh ) 
    
        
if __name__ == "__main__":
    main()
