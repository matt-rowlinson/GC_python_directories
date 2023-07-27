#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:30:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
import RowPy as rp
plt.style.use('seaborn-darkgrid')

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path+'/OutputDir/'):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', version='13.1.2'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    return ds

def site_data(ds, lat=16.9, lon=-24.9, z=0):
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1
    data = ds.isel( lon=x, lat=y, lev=z)
    return data

def find_model_output_for_site( site, rundirs, version, year='2017'):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        print( rundir )
        ds0 = get_data_as_xr(rundir, year=year, version=version)
        ds0 = site_data( ds0, lon=site['longitude'], lat=site['latitude'] )
        timeseries.append( ds0 )
    return timeseries

def plot(ds, GC_var, CV_var, ylabel):
    f, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(15,5),sharex=True)
    axes=[ax1,ax2,ax3,ax4,ax5,ax6]
    cs=['b','y','r']
    for n, ax in enumerate(axes):
        if CV_var[n]=='ALK4':
            iB = CV.get_from_merge(d['iso_butane']).fillna(0)
            nB = CV.get_from_merge(d['n_butane']).fillna(0)
            iP = CV.get_from_merge(d['iso_pentane']).fillna(0)
            nP = CV.get_from_merge(d['n_pentane']).fillna(0)
            df = iB + nB + iP + nP
        else:
            df = CV.get_from_merge(d[CV_var[n]])
        
        df = df.loc['2017'].resample( 'D' ).mean()
        ax.plot( df, c='k', label='Observations' )
        ax.plot( ds[0].time, ds[0][f'SpeciesConc_{GC_var[n]}'] * 1e12, color=cs[0], label='Base' )
        #if n < 2:
        #    ax.plot( ds[1].time, ds[1][f'SpeciesConc_{GC_var[n]}'] * 1e12, color=cs[1], label='Geo emissions' )
        ax.plot( ds[1].time, ds[1][f'SpeciesConc_{GC_var[n]}'] * 1e12, color=cs[2], label='Scaled CEDS' )

        ax.legend()
        ax.set_ylabel( f'{ylabel[n]} (pptv)' )
        if n > 3 :
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.tight_layout()
    return f 
##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    ## Set variables and read in GEOS-Chem data
    site    = sites['CVO']
    rundirs = ['geo_only','all_scaled']
    GC_var  = ['C2H6','C3H8','ALK4','BENZ','TOLU','PRPE']
    CV_var  = ['ethane','propane','ALK4','benzene','toluene','PRPE']
    ylabel  = ['Ethane','Propane','ALK4','Benzene','Toluene','Propene']
    ds = find_model_output_for_site(site, rundirs, version='13.4.0', year='2017')
    
    ## Make the plot and plot variable for each rundir
    f = plot( ds, GC_var, CV_var, ylabel) 
    plt.savefig( 'plots/TEST.figure_6.png', dpi=200 )
    plt.close()

if __name__ == "__main__":
    main()
