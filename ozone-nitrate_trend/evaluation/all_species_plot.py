#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=all_species_plot
#SBATCH --ntasks=1
#SBATCH --mem=36Gb
#SBATCH --partition=nodes
#SBATCH --time=10:00:00
#SBATCH --output=Logs/all_species_plots_%A.log

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
plt.style.use('seaborn-darkgrid')

def get_arguments():
    parser = argparse.ArgumentParser(description="Parse arguments to pass to GC processing scripts")
    parser.add_argument("-r", "--rundir", type=str, 
                        help='Name of desired GC rundir')
    parser.add_argument("-s", "--site", type=str,
                        default="CVO",
                        help="GAW site of interest")
    args=parser.parse_args()
    return args.rundir, args.site

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    print( len(file_list), "files to open")
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


##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():
    rundir, site = get_arguments()
    from CVAO_dict import CVAO_dict as d
    site = sites[site]
    ds = get_data_as_xr(f'nitrate_photol_control', year='201601')
    keys = sorted( ds.keys() )#[243:]
    print( len( keys ))

    #species_list = ['NH3']#NITS','NIT','NO','NO2','N','O3','SO2']
    #print( keys )
    sys.exit()
    
    for n, key in enumerate(keys):
        name = key.replace("SpeciesConc_","")
        #if name in species_list:
        if len(ds[key].shape)  < 4:
            continue
        print( key, n)
        data = ds[key]
        data = site_data( data, ds, lon=site['longitude'], lat=site['latitude'] ).resample(time='1MS').mean(dim='time')

        fig = plt.figure(n+1, figsize=(12,4))
        plt.plot( data['time'], data )
        plt.title(key)
        plt.ylabel( ds[key].units )

        plt.savefig( f'SpeciesConc_plots/{site["site_abbr"]}_{key}.png')
        plt.close()
        #else:
        #    pass

if __name__ == "__main__":
    main()
