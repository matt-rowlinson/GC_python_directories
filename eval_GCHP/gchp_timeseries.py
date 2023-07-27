#!/usr/bin/env python

# Imports
import gcpy.constants as gcon
import os
import numpy as np
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import xarray as xr
import glob
import warnings
import sys
sys.path.append('/users/mjr583/python_lib')
from CVAO_dict import CVAO_dict as d
import GC_tools as GC
import CVAO_tools as CV
# Tell matplotlib not to look for an X-window, as we are plotting to
# a file and not to the screen.  This will avoid some warning messages.
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

# Suppress harmless run-time warnings (mostly about underflow in division)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)


def find_files_in_dir(path, substrs):
    '''
    Returns a list of all files in a directory that match one or more
    substrings.

    Args:
    -----
        path : str
            Path to the directory in which to search for files.

        substrs : list of str
            List of substrings used in the search for files.

    Returns:
    --------
        file_list : list of str
            List of files in the directory (specified by path)
            that match all substrings (specified in substrs).
    '''
    
    # Initialize
    file_list = []

    # Walk through the given data directory.  Then for each file found,
    # add it to file_list if it matches text in search_list.
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))

    # Return an alphabetically sorted list of files
    file_list.sort()
    return file_list


def find_value_index(seq, val):
    '''
    Finds the index of a numpy array that is close to a value.

    Args:
    -----
        seq : numpy ndarray
            An array of numeric values.

        val : number
            The value to search for in seq.

    Returns:
    --------
        result : integer
            The index of seq that has a value closest to val.

    Remarks:
    --------
    This algorithm was found on this page:
    https://stackoverflow.com/questions/48900977/find-all-indexes-of-a-numpy-array-closest-to-a-value
    '''
    r = np.where(np.diff(np.sign(seq - val)) != 0)
    idx = r + (val - seq[r]) / (seq[r + np.ones_like(r)] - seq[r])
    idx = np.append(idx, np.where(seq == val))
    idx = np.sort(idx)
    result = np.round(idx)

    # NOTE: xarray needs integer values, so convert here!
    return int(result[0])


def read_geoschem_data(path, collections):
    '''
    Returns an xarray Dataset containing timeseries data.

    Args:
    -----
        path : str
            Directory path where GEOS-Chem diagnostic output
            files may be found.

        collections: list of str
            List of GEOS-Chem collections.  Files for these 
            collections will be read into the xarray Dataset.
            
    Returns:
    --------
        ds : xarray Dataset
            A Dataset object containing the GEOS-Chem diagnostic
            output corresponding to the collections that were
            specified.
    '''

    # Get a list of variables that GCPy should not read.
    # These are mostly variables introduced into GCHP with the MAPL v1.0.0
    # update.  These variables contain either repeated or non-standard
    # dimensions that can cause problems in xarray when combining datasets.
    skip_vars = gcon.skip_these_vars
    
    # Find all files in the given 
    file_list = find_files_in_dir(path, collections) 
    # Return a single xarray Dataset containing data from all files
    # NOTE: Need to add combine="nested" for xarray 0.15 and higher
    v = xr.__version__.split(".")
    #if int(v[0]) == 0 and int(v[1]) >= 25:
    #    return xr.open_mfdataset(file_list,
    #                             drop_variables=skip_vars,
    #                             combine="nested",
    #                             concat_dim=None)
    #else:
    print(file_list)
    return xr.open_mfdataset(file_list,
                                 drop_variables=skip_vars)


def plot_timeseries_data(ds, site_coords):
    '''
    Plots a timseries of data at a given (lat,lon) location.

    Args:
    -----
        ds : xarray Dataset
            Dataset containing GEOS-Chem timeseries data.

        site_coords : tuple
            Contains the coordinate (lat, lon) of a site location
            at which the timeseries data will be plotted.
    '''

    # ----------------------------------------------------------------------
    # Get the GEOS-Chem data for O3 and HNO3 corresponding to the
    # location of the observational station.  We will save these into
    # xarray DataArray objects, which we'll need for plotting.
    #
    # YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!
    # ----------------------------------------------------------------------
    
    # Find the indices corresponding to the site lon and lat
    #lat_idx = find_value_index(ds.lat.values, site_coords[0])
    #lon_idx = find_value_index(ds.lon.values, site_coords[1])
    variable='NO'
    # Save O3 from the first level (~60m height) (ppb) into a DataArray
    print('extract data')
    O3_L1 = ds['SpeciesConc_'+variable].isel(Xdim=12, Ydim=12,nf=5, lev=0)
    O3_L1 *= 1.0e12
    O3_L1.attrs['units'] = 'ppbv'
    
    a, lat, lon, lev, atime = GC.get_gc_var(rundir='tropchem_merra_4x5', variable=variable, version='12.9.3', year='2016')
    b, lat,lon,lev,btime = GC.get_gc_var(rundir='fullchem_4x5_LVOCfalse', variable=variable, version='GEOS-Chem',year='2016')

    y = find_value_index(lat, 16.9)
    x = find_value_index(lon, -24.9)

    a_time=[] ; b_time=[]
    for t in range(len(atime)):
        v=a[t,0,y,x]
        a_time.append(v)

    for t in range(len(btime)):
        v=b[t,0,y,x]
        b_time.append(v)
    import pandas as pd
    a=pd.DataFrame({'v12.9.3':a_time}, index=atime)
    b=pd.DataFrame({'v13.0.0':b_time}, index=btime)
    a=a['2016-01']
    b=b[a.index[0]:a.index[-1]]
    
    df=CV.get_from_merge(d[variable])
    df=df[a.index[0]:a.index[-1]]
    df=df.resample('D').mean()

    print('plotting')
    fig, ax = plt.subplots(1,1,figsize=(12,4))
    #O3_L1.plot.line(ax=ax, label='Ozone at Cape Verde (maybe)')
    ax.plot(a.index, O3_L1, label='GCHP')
    ax.plot(a.index,a, label='v12')
    ax.plot(b.index,b, label='v13')
    ax.plot(df.index, df.Value, 'k', label='CVAO')

    plt.ylabel(variable+' (ppbv)')
    plt.legend() 
    print('saving')
    plt.savefig('plots/CVO_lineplot_%s.png' %variable)

def main():
    '''
    Main program.
    '''


    # Path where the data files live
    # (YOU MUST EDIT THIS FOR YUR OWN PARTICULAR APPLICATION!)
    path_to_data = '/users/mjr583/scratch/GC/GCHP/rundirs/hold/'
    #ds = xr.open_dataset(path_to_data)

    # Get a list of files in the ConcAboveSfc and SpeciesConc collections
    # (YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!)
    collections = ['SpeciesConc']
    
    # Read GEOS-Chem data into an xarray Dataset
    ds = read_geoschem_data(path_to_data, collections)
    
    # Plot timeseries data at Centerville, AL (32.94N, 87.18W)
    # (YOU CAN EDIT THIS FOR YOUR OWN PARTICULAR APPLICATION!)
    site_coords = (16.9, -24.9)
    plot_timeseries_data(ds, site_coords)


if __name__ == "__main__":
    main()
