#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=global_ratios
#SBATCH --ntasks=1
#SBATCH --mem=8Gb
#SBATCH --partition=interactive
#SBATCH --time=00:45:00
#SBATCH --output=Logs/long_plots_%A.log
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

def get_data_as_xr(rundir, year=''):
    path=f'/mnt/lustre/users/mjr583/GC/13.1.2/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path, [f'SpeciesConc.{year}']), combine='by_coords' )
    return ds

def site_data(data, ds, lat=16.9, lon=-24.9, lev=0):
    x = rp.find_nearest(ds.lon, lon)
    y = rp.find_nearest(ds.lat, lat)
    if type(lev)==int:
        data = data.isel( lev=0, lon=x, lat=y )
    else:
        data = data.isel( lon=x, lat=y ).isel( lev=slice(lev)).mean(dim='lev', keep_attrs=True)
    return data

def load_observations( site, variable ):
    if site['site_name'] == 'CVAO':
        df = pd.read_csv( site['filepath'], index_col=0, dtype={"Airmass":str, "New_Airmass":str})
        #df.index = pd.to_datetime( df.index, format='%Y-%m-%d %H:%M:%S' )
        df = df[d[variable]['merge_name']]
    else:
        df = pd.read_csv( site['filepath'],index_col=0)
        df = df[variable]
    df.index = pd.to_datetime( df.index, format="%Y-%m-%d %H:%M:%S")
    return df


def loc_eval(x, b):
    loc_est = 0
    for i in enumerate(b): loc_est+=i[1]*(x**i[0])
    return(loc_est)


def loess(xvals, yvals, data, alpha, poly_degree=1):
    all_data = sorted(zip(data[xvals].tolist(), data[yvals].tolist()), key=lambda x: x[0])
    xvals, yvals = zip(*all_data)
    evalDF = pd.DataFrame(columns=['v','g'])
    n = len(xvals)
    m = n + 1
    q = int(np.floor(n * alpha) if alpha <= 1.0 else n)
    avg_interval = ((max(xvals)-min(xvals))/len(xvals))
    v_lb = min(xvals)-(.5*avg_interval)
    v_ub = (max(xvals)+(.5*avg_interval))
    v = enumerate(np.linspace(start=v_lb, stop=v_ub, num=m), start=1)
    xcols = [np.ones_like(xvals)]
    for j in range(1, (poly_degree + 1)):
        xcols.append([i ** j for i in xvals])
    X = np.vstack(xcols).T

    for i in v:
        #print(i)
        iterpos = i[0]
        iterval = i[1]
        iterdists = sorted([(j, np.abs(j-iterval)) for j in xvals], key=lambda x: x[1])
        _, raw_dists = zip(*iterdists)
        scale_fact = raw_dists[q-1]
        scaled_dists = [(j[0],(j[1]/scale_fact)) for j in iterdists]
        weights = [(j[0],((1-np.abs(j[1]**3))**3 if j[1]<=1 else 0)) for j in scaled_dists]
        _, weights      = zip(*sorted(weights,     key=lambda x: x[0]))
        _, raw_dists    = zip(*sorted(iterdists,   key=lambda x: x[0]))
        _, scaled_dists = zip(*sorted(scaled_dists,key=lambda x: x[0]))
        W         = np.diag(weights)
        b         = np.linalg.inv(X.T @ W @ X) @ (X.T @ W @ yvals)
        local_est = loc_eval(iterval, b)
        iterDF2   = pd.DataFrame({
            'v'  :[iterval],
            'g'  :[local_est]
            })
        evalDF = pd.concat([evalDF, iterDF2])
    evalDF = evalDF[['v','g']]
    return(evalDF)


def ignore_nans(df):
    idx=np.isfinite(df.Yvalue)
    Yvalue=df.Yvalue[idx]
    Xvalue=df.Xvalue[idx]
    index=df.index[idx]
    df=pd.DataFrame({'Xvalue':np.arange(len(Yvalue)), 'Yvalue':Yvalue }, index=index )
    return df


def linear_regression_model(df,X,Y, timestep='monthly'):
    ''' Guess of polynomial terms '''
    M = df.index.month

    z, p = np.polyfit(X, Y, 1)
    a = np.mean(Y.resample('A').mean())
    b = z
    g = 1.0
    d = 1.0
    M = 0
    def init_func(t,a,b,g,M,d):
        return a + b*t + g*np.cos(2*np.pi*M / 12) + d*np.sin(2*np.pi*M / 12)  

    guess = np.array([a, b, g, M, d])
    c,cov = curve_fit(init_func, X, Y, guess)
        
    def re_func(n,t,a,b,g,M,d,R):
        return a + b*t + g*np.cos(2*np.pi*M / 12) + d*np.sin(2*np.pi*M / 12)

    n = len(X)
    y = np.empty(n)
    R = np.empty(n)
    M = df.index.month

    for i in range(n):
        y[i] = re_func(i, X[i],c[0],c[1],c[2],M[i],c[4],R[i])
    var=y
    var = c[0] + c[1]*X #+ c[2]*X**2
    rmse = np.round(np.sqrt(mean_squared_error(Y,init_func( X, *c))),2)
    r2 = np.round(r2_score(Y,init_func( X, *c))*100,1)
    return y, var, z, rmse, r2, c, 


def get_anomaly(df):
    monmean=df.groupby(df.index.month).mean()
    anom = []
    for n in range(len(df)):
        nmonth=df.index[n].month
        anom.append( df[df.columns[1]][n] - monmean[monmean.columns[1]][nmonth] )
        
    anom = pd.DataFrame({'anomaly' : anom}, index=df.index)
    df=pd.DataFrame({'Xvalue':np.arange(len(anom.values)), 'Yvalue':anom.anomaly.values }, index=df.index )
    return df


def trend_statement(df):
    df = get_anomaly( df )
    output = linear_regression_model(df, df.Xvalue, df.Yvalue, timestep="M")
    trend_per_year=np.round(output[-1][1]*12,2)
    trend_per_decade=np.round(output[-1][1]*12*10, 2)
    print(f'O3 trend per year = {trend_per_year} ppb yr-1')
    return trend_per_year, output



##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    rundir, variables, version, site = get_arguments()
    print( rundir, variables, version, site)
    from CVAO_dict import CVAO_dict as d
    
    site = sites[site]
    variables=['O3']#,'CO','ethane'] 
    f, ax = plt.subplots(1,1,figsize=(10,4))
    for v in variables:
        data = pd.read_csv(site["save_name"]+'_O3.csv', index_col=0 )
        data.index= pd.to_datetime( data.index, format="%Y-%m-%d")

        df = load_observations(site, v)[:].resample("M").mean()
        print( site )
        if site['save_name']=='cvao':
            data=data['2006-10-01':]
            df['2009-07-01' : '2009-09-30'] = np.nan
            df = df['2006-10-01':'2017']
       
        ax.plot( data.index, data[site["save_name"]], c='g', label='v13.1.2', alpha=.5,zorder=2)
        data = pd.DataFrame({'Xvalue':np.arange(len(data[data.columns[0]].values)),'Yvalue' : data[data.columns[0]].values}, index=data.index)
        data = ignore_nans( data )
        evalDF = loess("Xvalue", "Yvalue", data = data, alpha=0.65, poly_degree=2)
        ax.plot(data.index, evalDF['g'][1:], color='g', linewidth= 2, zorder=12, label="_nolegend_")

        df.plot(c='k',zorder=1, alpha=.5, ax=ax)
        df.index.name = 'Xvalue'
        df=pd.DataFrame({'Xvalue':np.arange(len(df.values)), 'Yvalue':df.values }, index=df.index )
        df = ignore_nans(df)
        
        evalDF = loess("Xvalue", "Yvalue", data = df, alpha=0.65, poly_degree=2)
        ax.plot(df.index, evalDF['g'][1:], color='k', linewidth= 2, label="_nolegend_")

        gc_trend = trend_statement( data['1988':'2018'])
        obs_trend = trend_statement( df['1988':'2018'])

        plt.legend(loc=0)
        plt.ylabel(r"$O_3$ / ppb")
        plt.savefig( f"plots/{site['save_name']}_{v}.png" )
        plt.close()
        
if __name__ == "__main__":
    main()
