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
from sklearn.metrics import mean_squared_error

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
    return d, [x.Lon[-1].round(1), x.Lat[-1].round(1),int( x.Alt[-1].round(0)) ]

def process_datasets(ds, rundirs, year, d):#sites, year, var, gcvar):
    path='/mnt/lustre/users/mjr583/Ebas/data/'
    first=True
    for site in d['sites']:
        df = pd.read_csv( path+f'{site}_{d["var"]}.csv',index_col=0 )
        df.index = pd.to_datetime( df.index, format='%Y-%m-%d %H:%M')
        df = df.resample('D').mean()
        try:
            if "flag" in df.columns:
                df[df.flag != 0. ] = np.nan
            else:
                df[df['flag_'+d['gc']] != 0.] = np.nan
        except:
            pass

        d0, coords = model_data_for_comparison(ds[0], df, rundirs[0], year, d['gc'] )
        d1, coords = model_data_for_comparison(ds[1], df, rundirs[1], year, d['gc'] )
        
        df['site'] = site
        d0['site'] = site
        d1['site'] = site

        if first:
            first=False
            df_ = df.copy()
            d0_ = d0.copy()
            d1_ = d1.copy()
        else:
            df_  = pd.concat([df_,df])
            d0_  = pd.concat([d0_,d0])
            d1_  = pd.concat([d1_,d1])
    return df_, d0_, d1_

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


def plot_aromatic(ax,d,colours,df,d0,d1): 
    #### If LoD fixed to 10 ppt
    #df = df[df[d['form']] > 10. ]
    #d0 = d0[d0[d0.columns[0]] > 10.]
    #d1 = d1[d1[d1.columns[0]] > 10.]

    x0_=[] ; x1_=[] ; y0_=[] ; y1_=[]
    for site, marker in zip(d['sites'], d['markers']):
        x  = df[df.site == site]
        y0 = d0[d0.site == site]
        y1 = d1[d1.site == site]
        alt = int(np.nanmean(x.Alt.values))
        
        x  = x[d['form']].dropna() 
        y0 = y0.reindex( x.index )[y0.columns[0]]
        y1 = y1.reindex( x.index )[y1.columns[0]]

        ax.scatter( x, y0, c=colours[0], marker=marker, alpha=.4 )
        ax.scatter( x, y1, c=colours[1], marker=marker, alpha=.4 )
        ax.scatter( np.nan, 1, c='grey',marker=marker, label=f'{site} ({alt}m)' )

        x0 = x.reindex( y0.dropna().index )
        x1 = x.reindex( y1.dropna().index )
        x0_.append(x0) ; x1_.append(x1) ; y0_.append(y0) ; y1_.append(y1)
    
    x0 = pd.concat(x0_)
    x1 = pd.concat(x1_)
    y0 = pd.concat(y0_).dropna()
    y1 = pd.concat(y1_).dropna()
    
    x_fit, fit, a2 = get_odr( x0, y0 )
    nmb0 = 100 * np.sum( y0 - x0) / np.sum(x0)
    ax.plot( 10**x_fit, 10**fit, c='darkorange', zorder=200, label=f'Base (nmb={nmb0.round(2)}%)')
    
    x_fit, fit, a2 = get_odr( x1, y1 )
    nmb1 = 100 * np.sum( y1 - x1) / np.sum(x1)
    ax.plot( 10**x_fit, 10**fit, c='blue', zorder=200, label=f'Re-speciated (nmb={nmb1.round(2)}%)')
        
    if d['var']=='benzene':
        ax.legend(ncol=2, fontsize=6) 
    else:
        ax.legend(fontsize=7) 
    ax.set_xscale('log') ; ax.set_yscale('log') 
    ax.set_xlabel( f'Observed {d["var"].title()} (ppt)')
    ax.set_ylabel( f'Modelled {d["var"].title()} (ppt)')

    line  = matplotlib.lines.Line2D([0,1],[0,1], color='darkgrey', alpha=1., zorder=1)
    line.set_transform(ax.transAxes)
    ax.add_line(line)
    return
  

def var_dict():
    d = { 'benzene'     : {
              'var'     : 'benzene',
              'gc'      : 'BENZ',
              'form'    : 'C6H6',
              'sites'   : ['Zeppelin','Pallas','Hohenpeissenberg','Jungfraujoch','Zugspitze','Rigi','Cape_Verde'],
              'markers' : ['o','+','^','>','8','*','s'] },
          'toluene'   : {
              'var'   : 'toluene',
              'gc'    : 'TOLU',
              'form'  : 'C6H5CH3',
              'sites' : ['Zeppelin','Hohenpeissenberg','Rigi','Monte_Cimone','Cape_Verde'],
              'markers' : ['o','^','*','x','s'] },
          'xylene'   : {
              'var'   : 'xylene',
              'gc'    : 'XYLE',
              'form'  : 'XYLE',
              'sites' : ['Hohenpeissenberg','Rigi','Monte_Cimone'],
              'markers' : ['^','*','x'] }
              }
    return d
###----------------------------------------------------------------------------------------------------------------
def main():
    d = var_dict()
    rundirs = ['geo_2x25','all_2x25']
    years    = ['2016','2017']
    ds0 = get_data_as_xr(rundirs[0], years)
    ds1 = get_data_as_xr(rundirs[1], years)
    
    ## Make the plot
    fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(5,11))
    markers = ['o','+','^','>','8','*','s']
    colours = ['darkorange','blue']
    sub=str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

    ## Read, process and plot data for each aromatic 
    df, d0, d1 = process_datasets( [ds0,ds1], rundirs, years, d['benzene'])
    plot_aromatic(ax1, d['benzene'],colours, df,d0,d1)

    df, d0, d1 = process_datasets( [ds0,ds1], rundirs, years, d['toluene'])
    plot_aromatic(ax2, d['toluene'],colours, df,d0,d1)

    df, d0, d1 = process_datasets( [ds0,ds1], rundirs, years, d['xylene'])
    plot_aromatic(ax3, d['xylene'],colours, df,d0,d1)
    
    ## Set equal x/y-axis limits
    ax1.set_xlim( 5e-1,1e3 ) ; ax1.set_ylim( 5e-1,1e3 )
    ax2.set_xlim( 1e-2,1e3 ) ; ax2.set_ylim( 1e-2,1e3 )
    ax3.set_xlim( 5e-3,6e2 ) ; ax3.set_ylim( 5e-3,6e2 )
    
    ## LoD at 10 ppt
    #ax1.set_xlim( 7,1e3 ) ; ax1.set_ylim( 7,1e3 )
    #ax2.set_xlim( 7,1e3 ) ; ax2.set_ylim( 7,1e3 )
    #ax3.set_xlim( 7,6e2 ) ; ax3.set_ylim( 7,6e2 )




    plt.tight_layout()
    plt.savefig( f'plots/figure_Other-VOCs_GAW.png', dpi=300) 
    #plt.savefig( f'plots/figure_Other-VOCs_GAW_LoD-10ppt.png', dpi=300)    
    plt.close()
    
if __name__ == "__main__":
    main()
