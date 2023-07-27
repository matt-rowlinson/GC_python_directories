#!/usr/bin/env python3
import os
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def read_ems(rundir, variable, path='/mnt/lustre/users/mjr583/GC/14.1.0/rundirs/', year='2017'):
    ds = xr.open_mfdataset( find_file_list(f"{path}/{rundir}/", [f'HEMCO_diagnostics.{year}']), combine='by_coords' )
    ds = ds.mean('time')
    if 'lev' in ds.dims:
        ds = ds.sum(dim='lev')
    ds = ds[variable]
    return ds

##----------------------------------------MAIN SCRIPT----------------------------------------##
def main():

    ## Path to your HEMCO files and date string to read, can point to multiple files
    path = '/mnt/lustre/users/mjr583/GC/14.0.1/rundirs/'
    year='20170101'

    ## List the variables to be plotted
    variables = ['EmisNO_Anthro',
                 'EmisC2H6_Anthro',
                 'EmisC2H6_Total']
    
    ## Loop through the variables and create plot for each
    for v in variables:
        print( v )
        ## Change rundir name to your run directory
        base = read_ems( 'base_run',  v, path=path, year=year)
        dev  = read_ems( 'scale_all', v, path=path, year=year)
        
        ## Create the figure, plot each simulation and the difference ( dev-base)
        fig = plt.figure(figsize=(12,4))
        plots  = [base, dev, dev-base]
        labels = ['Base','Dev','Dev-Base']
        for n, ds in enumerate(plots):
            ax = plt.subplot(1,3,n+1, projection=ccrs.EqualEarth(central_longitude=0))
            if n==2:
                ds.plot.imshow( x='lon',y='lat', ax=ax, cmap='bwr', center=0.,  transform=ccrs.PlateCarree(),
                    cbar_kwargs={'orientation':'horizontal', 'label':labels[n]})
            else:
                ds.plot.imshow( x='lon',y='lat', ax=ax, transform=ccrs.PlateCarree(),
                    cbar_kwargs={'orientation':'horizontal', 'label':labels[n]})
            ax.coastlines()
        plt.tight_layout()
        plt.savefig(f'plots/hemco_diff_{v}.png')
        plt.close()

if __name__ == "__main__":
    main()
