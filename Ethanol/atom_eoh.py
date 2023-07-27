#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=cv-nitrate
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
#SBATCH --partition=interactive
#SBATCH --time=00:45:00
#SBATCH --output=Logs/cv-nitrate_%A.log
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
import glob
plt.style.use('seaborn-darkgrid')
sys.path.append('/users/mjr583/cvao')
import numpy as np
import seaborn as sns
import datetime

var_dict = { 'NO3' : {
                'mm' : 62.0049 },
             'NH4' : {
                'mm' : 18.0400 },
             'SO4' : {
                'mm' : 96.06 },
             'Na' : {
                'mm' : 22.989769 },
             'Ca2': {
                'mm' : 40.00 },
             'EOH' : {
                'mm' : 46.06844 }
        }

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    print( file_list )
    return file_list

def get_data_as_xr(rundir, year='', variable=False, version='13.1.2'):
    path=f'/mnt/lustre/users/mjr583/GC/{version}/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    return ds

def site_data(ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)
    y = rp.find_nearest(ds.lat, lat)
    if type(lev)==int:
        data = ds.isel( lon=x, lat=y, lev=lev )
    else:
        data = ds.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

flight_dict = {
        'ATom1' : { 'date':'2016-08', 'modeldate':'201708' },
        'ATom2' : { 'date':'2017-02', 'modeldate':'201702' },
        'ATom3' : { 'date':'2017-10', 'modeldate':'201710' },
        'ATom4' : { 'date':'2018-05', 'modeldate':'201705' },          
            }

def read_gc(rundir, label, flight):
    ds = get_data_as_xr(rundir, year=flight_dict[flight]['modeldate'])
    # Read GEOS-Chem data
    ds0 = ds['SpeciesConc_EOH'] * 1e6 * 10
    ds0 = ds0.isel( lat=slice(25,31) ).mean( dim='lat' )
    ds0 = ds0.isel( lon=slice(30,34) ).mean( dim='lon' )

    alt=pd.read_csv('/users/mjr583/GC/info_files/GC_72_vertical_levels.csv')['Altitude (km)'] * 1e3  
    aa=[] ; alts=[]
    for lev in range(30):
        a = ds0[:,lev] #var_dict['EOH']['mm'] * ( ds0[:,lev] / ( 22.41 * ( t[:,lev] / 273 ) * ( 1013 / p[:,lev] ) ) )
        a_alt = np.repeat( alt[lev], len(ds0) )
        alts.append( a_alt )#.values )
        aa.append( a.values )
    
    aa = np.array(aa).flatten()
    alts = np.array(alts).flatten()
    ds = pd.DataFrame( {'Ethanol_GMI':aa, 'G_ALT':alts } )
    ds['Model']=label
    return ds

def read_atom(flight):
    infile = '/mnt/lustre/groups/chem-acm-2018/shared_data/NASA/ATom/v2020-07-05/Mor.all.at1234.2020-07-05.tbl'
    df = pd.read_table(infile, comment='#', delimiter=',')
    days= pd.to_datetime( df['YYYYMMDD'] , format="%Y%m%d" )
    times= df['UTC_Mean'] 
    timex=[] 
    for n, i in enumerate(days):
        timex.append( i + datetime.timedelta( 0, times[n] ) )
    df.index = pd.to_datetime( timex, format='%YYYY-%m-%d %H:%M:%S')
    df = df[df['Ethanol_GMI'] > 0.]
    df = df.loc[flight_dict[flight]['date']]
    df = df[['G_LAT','G_LONG','G_ALT','Ethanol_GMI']]
    df = df.loc[(df.G_LAT > 10.) & (df.G_LAT < 25.) ]
    df = df.loc[(df.G_LONG > -30.) & (df.G_LONG < -15.) ]
    df['Model'] = flight
    return df


def main():
    flights=['ATom1','ATom2','ATom3','ATom4']
    for flight in flights:
        print( flight )
        # Get correct T and p
        path='/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/nitrate_photol_control_DUST/'
        ds = xr.open_mfdataset( find_file_list(path, [f'StateMet.{flight_dict[flight]["modeldate"]}']), combine='by_coords' )
        t = ds['Met_T'].isel( lat=slice(25,31) ).mean( dim='lat' ).isel( lon=slice(30,34) ).mean( dim='lon' )
        p = ds['Met_PMIDDRY'].isel( lat=slice(25,31) ).mean( dim='lat' ).isel( lon=slice(30,34) ).mean( dim='lon' )
         
        ds0 = read_gc('ceds_only', 'Base CEDS', flight)
        ds1 = read_gc('new_scale_all_vocs', 'Scaled CEDS', flight)
        ds2 = read_gc('asia-emep_scale_all', 'Upper EOH scale', flight)
        
        # Read ATom data
        df = read_atom(flight)
        '''
        infile = '/mnt/lustre/groups/chem-acm-2018/shared_data/NASA/ATom/v2020-07-05/Mor.all.at1234.2020-07-05.tbl'
        df = pd.read_table(infile, comment='#', delimiter=',')
        days= pd.to_datetime( df['YYYYMMDD'] , format="%Y%m%d" )
        times= df['UTC_Mean'] 
        timex=[] 
        for n, i in enumerate(days):
            timex.append( i + datetime.timedelta( 0, times[n] ) )
        df.index = pd.to_datetime( timex, format='%YYYY-%m-%d %H:%M:%S')
        df = df[df['Ethanol_GMI'] > 0.]
        df = df.loc[flight_dict[flight]['date']]
        df = df[['G_LAT','G_LONG','G_ALT','Ethanol_GMI']]
        df = df.loc[(df.G_LAT > 10.) & (df.G_LAT < 25.) ]
        df = df.loc[(df.G_LONG > -30.) & (df.G_LONG < -15.) ]
        df['Model'] = flight
        '''
        
        # Merge the two data sources into one dataframe 
        x = df['G_ALT'].append( ds0.G_ALT ).append( ds1.G_ALT ).append( ds2.G_ALT )
        y = df['Model'].append( ds0.Model ).append( ds1.Model ).append( ds2.Model )
        z = df['Ethanol_GMI'].append( ds0['Ethanol_GMI'] ).append( ds1['Ethanol_GMI'] ).append( ds2['Ethanol_GMI'] )
        ds = pd.DataFrame( {'ALTITUDE': x.values, 'Data':y.values,'EOH':z.values } )

        bins = [0, 1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3]#, 1e4]
        ds['binned'] = np.searchsorted(bins, ds['ALTITUDE'].values) 

        # PLOTTING 
        f,ax=plt.subplots()
        sns.boxplot(x="EOH", y="binned", width=.5, hue="Data",palette=['darkgrey',"c", "b", "r"], 
                orient='h',data=ds, ax=ax, showfliers=False)# labels=["Metric", "Length"])
        
        ax.invert_yaxis()
        ax.set_ylabel('Altitude (km)')
        ax.set_xlabel('EOH ppb')
        plt.tight_layout()
        plt.savefig( f'plots/EOH_ATom_{flight}.png' )
        plt.close()

if __name__=="__main__":
    main()

