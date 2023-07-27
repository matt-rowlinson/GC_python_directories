#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=2Gb
#SBATCH --partition=interactive
#SBATCH --time=00:05:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

import glob
import xarray as xr
import os 
import numpy as np
from matplotlib.offsetbox import AnchoredText
import matplotlib.lines as mlines
import matplotlib.dates as mdates

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, years, lev=0, version='14.0.1', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    substrs = ['SpeciesConc.'+s for s in years]
    ds = xr.open_mfdataset( find_file_list(path, substrs), combine='by_coords' )
    return ds

def find_nearest_3d(ds,lon_value, lat_value, alt_value):
    gc_alts = pd.read_csv( '/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)'] * 1e3
    x_idx=(np.abs( ds.lon - float(lon_value) )).argmin()
    y_idx=(np.abs( ds.lat - float(lat_value) )).argmin()
    z_idx=(np.abs( gc_alts.values - float(alt_value) )).argmin()
    return x_idx, y_idx, z_idx

def model_data_for_comparison(ds, df, rundir, year, var ):
    #ds = get_data_as_xr(rundir, year)
    d = ds[f'SpeciesConc_{var}']
    
    x = df.dropna()
    x_idx, y_idx, z_idx = find_nearest_3d( ds, x.Lon[-1], x.Lat[-1], x.Alt[-1] )
    d = d.isel( lon=x_idx, lat=y_idx, lev=z_idx )
    d = pd.DataFrame( {f'{rundir}': d.values * 1e12}, index=d.time ).resample('D').mean()
    return d, [x.Lon[-1].round(1), x.Lat[-1].round(1), x.Alt[-1].round(1) ]

def process_datasets(ds, rundirs, site, year, var, gcvar):
    path='/mnt/lustre/users/mjr583/Ebas/data/'
    df = pd.read_csv( path+f'{site}_{var}.csv',index_col=0 )
    df.index = pd.to_datetime( df.index, format='%Y-%m-%d %H:%M')
    df = df.resample('D').mean()
    if "flag" in df.columns:
        df[df.flag != 0. ] = np.nan
    else:
        df[df['flag_'+gcvar] != 0.] = np.nan
    
    d0, coords = model_data_for_comparison(ds[0], df, rundirs[0], year, gcvar ) 
    d1, coords = model_data_for_comparison(ds[1], df, rundirs[1], year, gcvar )
    return df, d0, d1,coords

def func(p, x):
    a, b = p
    return a * x + b

def get_odr(x,y):
    from scipy import odr
    quad_model = odr.Model(func) 

    data = odr.RealData( np.log10(x),np.log10(y) )
    odr  = odr.ODR( data, quad_model, beta0=[0., 1.,] ) 
    out  = odr.run()
    popt = out.beta
    perr = out.sd_beta
    #print( popt )

    for  i in range(len(popt)):
        nstd=10
        popt_up = popt + nstd * perr
        popt_dw = popt - nstd * perr
        x_fit = np.linspace( np.log10(x).min(), np.log10(x).max(), len(x) )
        fit = func( popt, x_fit )
        return x_fit, fit, popt[0] 

###----------------------------------------------------------------------------------------------------------------
def main():
    # Ethane sites
    rundirs = ['geo_2x25','all_2x25']
    years    = ['2016','2017']
    var     = 'benzene'
    gcvar   = 'BENZ'
    form    = 'C6H6'
    sites = ['Zeppelin','Pallas','Hohenpeissenberg','Jungfraujoch','Zugspitze','Rigi','Monte_Cimone','Cape_Verde']
    ds0 = get_data_as_xr(rundirs[0], years)
    ds1 = get_data_as_xr(rundirs[1], years)
    
    fig = plt.figure(figsize=(11,15.2), tight_layout=True)
    gs  = matplotlib.gridspec.GridSpec(8, 4)
    
    ax1 = fig.add_subplot(gs[0, :-1])
    ax2 = fig.add_subplot(gs[1, :-1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2, :-1], sharex=ax1)
    ax4 = fig.add_subplot(gs[3, :-1], sharex=ax1)
    ax5 = fig.add_subplot(gs[4, :-1], sharex=ax1)
    ax6 = fig.add_subplot(gs[5, :-1], sharex=ax1)
    ax7 = fig.add_subplot(gs[6, :-1], sharex=ax1)
    ax8 = fig.add_subplot(gs[7, :-1], sharex=ax1)

    ax9  = fig.add_subplot(  gs[0, 3])
    ax10 = fig.add_subplot(  gs[1, 3])
    ax11 = fig.add_subplot(  gs[2, 3])
    ax12 = fig.add_subplot(  gs[3, 3])
    ax13 = fig.add_subplot(  gs[4, 3])
    ax14 = fig.add_subplot(  gs[5, 3])
    ax15 = fig.add_subplot(  gs[6, 3])
    ax16 = fig.add_subplot(  gs[7, 3])

    axes  = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
    axes2 = [ax9,ax10,ax11,ax12,ax13,ax14,ax15,ax16]

    colours = ['blue','darkorange'][::-1]
    sub=str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    Max = 900
    Min = 0.5
    for site, ax, ax_ in zip(sites, axes, axes2):
        print( site )
        df, d0, d1, coord = process_datasets( [ds0,ds1], rundirs, site, years, var, gcvar)
        ax.plot( d0[rundirs[0]], color=colours[0], label='Base' )
        ax.plot( d1[rundirs[1]], color=colours[1], label='Scaled CEDS emissions' )
        print( d1[rundirs[1]].min())
        print( d0[rundirs[0]].min())
        ax.plot( df[form], c='k', label='Observations' )

        textstr = (f'{site.replace("_"," ")}\n{coord[0]}$^\circ$W, {coord[1]}$^\circ$N\n{coord[2]}m') 
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5, pad=0.2)
        ax.text(0.05, 0.05, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='bottom', bbox=props)

        ax.set_ylabel( f'{gcvar.translate(sub)} (pptv)' )
        ax.set_ylim(bottom=Min, top=Max)
        ax.set_yscale('log')

        line  = matplotlib.lines.Line2D([0,1],[0,1], color='darkgrey', alpha=1., zorder=10)
        line.set_transform(ax_.transAxes)
        ax_.add_line(line)
        
        y = d0[rundirs[0]].reindex( df.dropna().index )
        ax_.scatter( df[form].dropna(),y, s=4, c='darkorange', alpha=.2, zorder=100)
        x_fit, fit, a2 = get_odr( df[form].dropna(), y )
        ax_.plot( 10**x_fit, 10**fit, c='darkorange', zorder=200)

        y = d1[rundirs[1]].reindex( df.dropna().index )
        ax_.scatter( df[form].dropna(),y, s=4, c='blue', alpha=.2, zorder=200)
        x_fit, fit, a2 = get_odr( df[form].dropna(), y )
        ax_.plot( 10**x_fit, 10**fit, c='b', zorder=200)

        ax_.set_xlim(left=Min, right=Max) ; ax_.set_ylim(bottom=Min, top=Max)
        ax_.set_xscale('log') ; ax_.set_yscale('log')
        
        ax_.set_ylabel( f'Model {gcvar.translate(sub)} (pptv)' )

        if site != sites[-1]:
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax_.get_xticklabels(), visible=False)

    ax1.legend(loc="upper center", bbox_to_anchor=(0.5,1.2), ncol=3, fontsize=12)
    ax16.set_xlabel( f'Observed {gcvar.translate(sub)} (pptv)' )
    myFmt = matplotlib.dates.DateFormatter('%m-%Y')
    ax.xaxis.set_major_formatter(myFmt)

    plt.tight_layout()
    plt.savefig( f'plots/figure_GAW_sites_{var}.png', dpi=300)   
    plt.close()

if __name__ == "__main__":
    main()
