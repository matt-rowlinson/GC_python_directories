#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=jscale_time
#SBATCH --ntasks=1
#SBATCH --mem=5Gb
#SBATCH --partition=nodes
#SBATCH --time=01:30:00
#SBATCH --output=Logs/jscale_time_plots_%A.log
import sys
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
from CVAO_dict import CVAO_dict as d
import RowPy as rp
plt.style.use('seaborn-darkgrid')
import os
import xarray as xr
import matplotlib.lines as mlines
from scipy import odr
import matplotlib.lines as mlines
plt.style.use('seaborn-darkgrid')

def func(p, x):

    a, b, = p
    return a * x + b

# Model object
quad_model = odr.Model(func)

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path+'/OutputDir/'):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, year='', collection='SpeciesConc', variable=False, Chem=False, version='13.1.2'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    if variable=='OH' or variable=='HO2':
        collection='ConcAfterChem'
        Chem=True
    try:
        if year=='201001':
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
        else:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )

    except:
        try:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}']), combine='by_coords' )
        except:
            ds = xr.open_mfdataset( find_file_list(path, [f'{collection}.{year}'])[:-1], combine='by_coords' )

    if variable:
        if variable=='NO3':
            try:
                ds = ( ds[f'SpeciesConc_NIT'] + ds[f'SpeciesConc_NITs'] + ds[f'SpeciesConc_NITD1'] 
                     + ds[f'SpeciesConc_NITD2']+ ds[f'SpeciesConc_NITD3'] + ds[f'SpeciesConc_NITD4']   ).isel(lev=0) * 1e9
            except:
                ds = ( ds[f'SpeciesConc_NIT'] + ds[f'SpeciesConc_NITs'] ).isel(lev=0) * 1e9
            ds = ds * ( 62.0049 / 24.45 )
    return ds

def site_data(ds1, ds, lat=16.9, lon=-24.9, lev=0):
    #if site['save_name']=='cvao':
    #    x = rp.find_nearest(ds.lon, lon)-1
    #    y = rp.find_nearest(ds.lat, lat)+1
    #else:
    x = rp.find_nearest(ds.lon, lon)-1
    y = rp.find_nearest(ds.lat, lat)+1

    if type(lev)==int:
        data = ds1.isel( lon=x, lat=y )
    else:
        data = ds1.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def find_model_output_for_site(ds, rundirs, versions, v, years='2017'):
    timeseries=[] ; add_lowes=[]
    for n, rundir in enumerate(rundirs):
        #print( rundir )
        if len( versions ) == 1:
            version=versions[0]
        else:
            version=versions[n]
        ds0 = get_data_as_xr(rundir, year=years, variable=v, version=versions[0])

        ds0 = site_data( ds0, ds )
        ds0 = pd.DataFrame({'CVAO':ds0.values}, index=ds0['time'].values)#.resample('M').mean()
        #if site['save_name']=='cvao':
        #    ds0=ds0['2006-10-01':]
        timeseries.append( ds0 )
    return timeseries

def main():
    #global d, site, rundirs, v
    sys.path.append('/users/mjr583/cvao/')
    from lowess_function import loess
    var='NO3'
    alpha=.35 ; poly=1
    
    df0 = pd.read_csv('/mnt/lustre/users/mjr583/NCAS_CVAO/CVAO_datasets/Nitrate_Ammonium_CVAO12-19.csv',
                    index_col=1) 

    
    df0.index = pd.to_datetime( df0.index, format="%d/%m/%Y %H:%M")
    df0 = df0.dropna()
    df0=df0[df0.columns[2]]
    df0=pd.DataFrame(df0)
    df0.columns=['Value']
    df0=df0['2018']

    
    df2 = pd.DataFrame( {"XValue" : np.arange(len(df0.Value.values)), "YValue" : df0.Value.values   }, index=df0.index  )
    df2 = df2.resample('W').mean()
    df2=df2.dropna()
    eval_DF = loess("XValue", 'YValue', data = df2, alpha=alpha, poly_degree=poly )

    ds = get_data_as_xr('nitrate_photol_control', year='201001')
    dsy = find_model_output_for_site(ds, ['dev_new_base'], ['13.1.2'], var, years='2018')[0]['2012':'2019']
    ds1 = pd.DataFrame( {"XValue" : np.arange(len(dsy.values)), "YValue" : dsy.CVAO.values   }, index=dsy.index  )
    ds1 = ds1.resample('W').mean()
    eval_DF2 = loess("XValue", 'YValue', data = ds1, alpha=alpha, poly_degree=poly )

    dsx = find_model_output_for_site(ds, ['dev_both'], ['13.1.2'], var, years='2018')[0]['2012':'2019']
    ds2 = pd.DataFrame( {"XValue" : np.arange(len(dsx.values)), "YValue" : dsx.CVAO.values   }, index=dsx.index  )
    ds2 = ds2.resample('W').mean()
    eval_DF3 = loess("XValue", 'YValue', data = ds2, alpha=alpha, poly_degree=poly )
    
    df0.index = pd.to_datetime(df0.index.strftime('%Y-%m-%d'), format='%Y-%m-%d')
    dsy.index = pd.to_datetime(dsy.index.strftime('%Y-%m-%d'), format='%Y-%m-%d')
    dsx.index = pd.to_datetime(dsy.index.strftime('%Y-%m-%d'), format='%Y-%m-%d')
    dsy = dsy.reindex( df0.index )
    dsx = dsx.reindex( df0.index )
    
    x=df0
    y=dsy
    z=dsx
    x0=x
    y0=y
    z0=z
    from scipy import odr
    # Model object
    quad_model = odr.Model(func)

    # Create a RealData object
    data = odr.RealData(np.log10(x), np.log10(y))
    # Set up ODR with the model and data.
    odr = odr.ODR(data, quad_model, beta0=[0., 1.])
    # Run the regression.
    out = odr.run()

    # Create a RealData object
    data2 = odr.RealData(np.log10(x), np.log10(z))
    # Set up ODR with the model and data.
    odr2 = odr.ODR(data, quad_model, beta0=[0., 1.])
    # Run the regression.
    out2 = odr.run()


    #print fit parameters and 1-sigma estimates
    popt = out.beta
    perr = out.sd_beta
    #print('fit parameter 1-sigma error')
    #print('———————————–')
    for i in range(len(popt)):
        print(str(popt[i])+' +- '+str(perr[i]))

        # prepare confidence level curves
        nstd = 10. # to draw 5-sigma intervals
        popt_up = popt + nstd * perr
        popt_dw = popt - nstd * perr

        x_fit = np.linspace(min(np.log10(x)), max(np.log10(x)), len(x))
        fit = func(popt, x_fit)
        fit_up = func(popt_up, x_fit)
        fit_dw= func(popt_dw, x_fit)
    
    colors=['#bae4bc','#7bccc4','#2b8cbe']
    f, ax = plt.subplots( 1,1, figsize=(6,6))

    line = mlines.Line2D([0, 1], [0, 1], color='darkgrey', alpha=1.)
    transform = ax.transAxes
    line.set_transform(transform)
    ax.add_line(line)
    
    ax.scatter( df0.values, dsy.values, c='g', label='Base' )
    ax.scatter( df0.values, dsx.values, c='orange', label='NIThv')

    #errorbar(x, y, yerr=y, xerr=x, hold=True, ecolor='k', fmt='none', label='data')
    #ax.plot(10**x_fit, 10**fit, c=colors[2], lw=2, label='best fit curve')

    ax.set_ylabel( 'Model $NO_3^-$ (ug $m^{-3}$)')
    ax.set_xlabel( 'Observed $NO_3^-$ (ug $m^{-3}$)')
    ax.set_xlim(left=0., right=2.6)
    ax.set_ylim(bottom=0., top=2.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig( 'plots/CVAO_scatter.png' )
    plt.close()
    
    
    #print( x.max(), x.min())
    #print( y.max(), y.min())
    ##errorbar(x, y, yerr=y, xerr=x, hold=True, ecolor='k', fmt='none', label='data')
    #ax.plot(10**x_fit, 10**fit, c=colors[2], lw=2, label='best fit curve')
    #print( 10**x_fit.max(), 10**x_fit.min())
    #print('!!!',  np.exp(x_fit).max(), np.exp(x_fit).min())

    #plot(x0, y0, 'k-', lw=2, label='True curve')
    #ax.fill_between(x_fit, fit_up, fit_dw, alpha=.25, label='5-sigma interval')

    #Max = np.max( [ np.nanmax( y ), np.nanmax( x ) ] )
    #Pc  = np.max( [ np.nanpercentile( y, 97 ), np.nanpercentile( x, 97 ) ] )

    #ax.set_ylabel( 'Model $NO_3^-$ (ug $m^{-3}$)')
    #ax.set_xlabel( 'Observed $NO_3^-$ (ug $m^{-3}$)')
    #ax.set_xlim(left=8e-4, right=2e1)
    #ax.set_ylim(bottom=8e-4, top=2e1)
    #ax.set_yscale('log')
    #ax.set_xscale('log')

    #plt.tight_layout()
    #plt.savefig( 'plots/CVAO_scatter.png' )
    #plt.close()


    '''

    ds25 = df0.resample('M').quantile(.25)
    ds75 = df0.resample('M').quantile(.75)
    ax1.fill_between(ds25.index, ds25.Value.values, ds75.Value.values, color='k', alpha=.2)

    #ax1.scatter(df2.index, df2["YValue"], color="k", marker="o", alpha=0.3,s=5, label="_nolegend_")
    #ax1.plot(df0.Value.resample('M').mean().index, df0.Value.resample('M').mean(), color="k", alpha=0.7, label="_nolegend_")
    #ax1.plot(df2.index, eval_DF['g'][1:], color='k', linewidth= 3, label="CVAO")

    #ax1.scatter( ds.index, ds.CVAO, color="green", marker="o", alpha=0.3, s=5, label="_nolegend_")
    #ax1.plot(ds.resample('M').mean().index, ds.resample('M').mean(), color="green", alpha=0.7, label="_nolegend_")
    ds25 = dsy.resample('M').quantile(.25)
    ds75 = dsy.resample('M').quantile(.75)
    ax1.fill_between(ds25.index, ds25.CVAO.values, ds75.CVAO.values, color='green', alpha=.2)

    ds25 = dsx.resample('M').quantile(.25)
    ds75 = dsx.resample('M').quantile(.75)
    ax1.fill_between(ds25.index, ds25.CVAO.values, ds75.CVAO.values, color='orange', alpha=.2)
    
    
    ax1.plot(df0.Value.resample('M').mean().index, df0.Value.resample('M').mean(), color="k", alpha=0.7, label="CVAO")
    ax1.plot(dsy.resample('M').mean().index, dsy.resample('M').mean(), color="green", alpha=0.7, label="GC_Base")
    ax1.plot(dsx.resample('M').mean().index, dsx.resample('M').mean(), color="orange", alpha=0.7, label="GC_NIThv")


    #ax1.plot(df2.index, eval_DF['g'][1:], color='k', linewidth= 2.5, label="CVAO")
    #ax1.plot(ds1.index, eval_DF2['g'][1:], color='green', linewidth= 2.5, label="GEOS-Chem")
    #ax1.plot(ds2.index, eval_DF3['g'][1:], color='orange', linewidth= 2.5, label="GEOS-Chem_NIThv")



    plt.xlabel(None), plt.ylabel(f"{var} / ug m-3")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'plots/cvao_{var}_lowess.png')
    plt.close()  
    '''

        
if __name__ == "__main__":
    main()
