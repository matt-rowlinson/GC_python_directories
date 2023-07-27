#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NOxtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=nodes
#SBATCH --time=00:45:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
from sites_dicts import GAW_dict as sites
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
plt.style.use('seaborn-darkgrid')
import os
import xarray as xr

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path+'/OutputDir/'):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def site_data(ds, site):
    x = rp.find_nearest(ds['lon'], site['longitude'])
    y = rp.find_nearest(ds['lat'], site['latitude'])
    if site['site_name']=='CVAO':
        x= x+1
        y= y-1
    data = ds.isel( lon=x, lat=y, lev=0 )
    return data

def get_data_as_xr(site, rundir, years, variable, version='13.1.2'):
    print(rundir)
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    if type(years)==list:
        list_of_lists=[]
        for year in years:
            path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}{year}'
            y1 = str( int(year) -1 )
            #print( path, 'SpeciesConc.{y1}' )
            #sys.exit()
            files = find_file_list(path, [f'SpeciesConc.{y1}'])
            list_of_lists.append( files )
        files = [x for xs in list_of_lists for x in xs]
        #print( files ) 
    first=True
    for f in files:
        try:
            d = xr.open_dataset( f )
        except: 
            continue
        
        ddd = site_data(d, site)
        ddd = ddd[f'SpeciesConc_{variable}']*1e9

        #dd = d['SpeciesConc_O3'].isel(lev=0) * 1e9
        #ddd = dd.isel( lon=31 ).isel( lat=27 )
        
        if first:
            df = pd.DataFrame( {'O3': ddd.values}, index=d['time'] )
            first=False
        else:
            af = pd.DataFrame( {'O3': ddd.values}, index=d['time'] ) 
            df = pd.concat([df, af])
        
    return df 

def get_model_as_dataframe(rundirs, variable, version="13.1.2", collection="SpeciesConc", year="2017", x=-24.9, y=16.9, z=0):  
    dfs=[]
    for a in rundirs:
        print(a)
        var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, collection=collection, version='13.1.2', year=year)
        a_time=[] 
        idy = rp.find_nearest(lat, y)
        idx = rp.find_nearest(lon, x)
        for t in range(len(atime)):
            v=var[t,z,idy,idx]
            a_time.append(v)
        a=pd.DataFrame({a:a_time}, index=atime)
        dfs.append(a)
    df=pd.concat(dfs, axis=1, ignore_index=False)
    return df

def trendline(df):
    XX=np.arange(len(df))
    idx=np.isfinite(df)
    Y=df[idx]
    X=XX[idx]
    time=df.index[idx]

    z = np.polyfit(X, Y, 1)
    p = np.poly1d(z)
    z = z[0] * 12. * 10.
    return idx, z, p, X

def dataframe_by_year(year):
    rundirs=[f'nitrate_photol_control']
    try:
        df = get_model_as_dataframe( rundirs, variable="O3", year=year )
    except:
        df = get_model_as_dataframe( rundirs, variable="O3", year='' )
    df.columns=['No_dust']#,'Base','j100']
    return df

def plot_by_scale(ax, df, label, c='r', first=False):
    idx, z, p, X = trendline( df )
    df=df.dropna()
    print(df)
    if first:
        ax.scatter( df.index, df, c=c,label=f'{label}  (y=%.2f ppb decade$^{-1}$)' %z)
        ax.plot( df.index, df, '--', c=c )
    else:
        ax.scatter( df.index, df, c=c,label=f'{label}  (y=%.2f)' %z)
        ax.plot( df.index, df, '--', c=c )

    #ax.plot( df.index[idx], p(X), '--', c=c)
    return

def trend(ax, df, label, c='r', first=False):
    idx, z, p, X = trendline( df )
    ax.plot( df.index[idx], p(X), '--', c=c)
    return

def load_observations( site, variable ):
    if site['site_name'] == 'CVAO':
        df = pd.read_csv( site['filepath'], index_col=0, dtype={"Airmass":str, "New_Airmass":str}, low_memory=False)
        if variable=='NOx':
            df = df['NO_pptV'] + df['NO2_pptV'] 
        else:
            df = df[d[variable]['merge_name']]
    elif site['site_name'] == 'Cape Grim':
        df = pd.read_csv( site['filepath']+site[f'{variable}_file'],index_col=0)
        df = df['Value']
        df = df / 1.9957
        df.index = pd.to_datetime( df.index, format="%d/%m/%Y %H:%M")
    else:
        df = pd.read_csv( site['filepath']+site[f'{variable}_file'],index_col=0)
        df = df['Value']
    if site['unit_conv'] == True:
        df = df / 1.96
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df


def main():
    site = sites['CVO']

    base = get_data_as_xr( site, 'base_run_', years=['1980','1990','2000','2005','2010','2015'],variable='O3' )
    j100 = get_data_as_xr( site, 'j100_run_', years=['1980','1990','2000','2005','2010','2015'],variable='O3' )

    m = pd.concat([base, j100], axis=1).resample( 'Y' ).mean()
    m.columns = ['Base','j100']
    print( m )
   
    df_1980 = dataframe_by_year('1980')
    df_2000 = dataframe_by_year('2000')
    df_2005 = dataframe_by_year('2005')
    df_2010 = dataframe_by_year('2010')
    df_2017 = dataframe_by_year('2017')

    merged = pd.concat( [df_1980, df_2000, df_2005, df_2010, df_2017 ] )
    merged=merged.resample('Y').mean() * 1e9
    

    f,ax=plt.subplots(figsize=(14,5))
    colors=['#1b9e77','#d95f02','#7570b3','#e7298a','k']
    #plot_by_scale(ax, merged['No_dust'], label='No Dust', c=colors[0], first=True)
    plot_by_scale(ax, m['Base'],    label='Base',  c=colors[1])
    plot_by_scale(ax, m['j100'],    label=r'J$_{100scale}$', c=colors[3])

    df = load_observations(site, 'O3')[:].resample("Y").mean()
    plot_by_scale(ax, df, label='CVAO', c='k')

    split_trend=False
    if split_trend:
        #trend(ax, m['Base'][:'2000'],    label='__Base',  c=colors[2])
        #trend(ax, m['Base']['2000':],    label='__Base',  c=colors[2])
        #trend(ax, m['j100'][:'2000'],    label='__Base',  c=colors[4])
        #trend(ax, m['j100']['2000':],    label='__Base',  c=colors[4])
        ax.set_ylabel( r'$O_3$ (ppbv)')    
        plt.legend()
        plt.tight_layout()
        plt.savefig('plots/trends_by_Jscaling_split.png')
        plt.close()
    else:
        ax.set_ylabel( r'$O_3$ (ppbv)')    
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'plots/ANNMEANS_{site["save_name"]}_trends_by_Jscaling.png')
        plt.close()



if __name__=="__main__":
    main()
