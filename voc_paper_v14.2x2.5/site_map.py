#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=site_map
#SBATCH --ntasks=1
#SBATCH --mem=500Mb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=sitemap_%A.log
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import cartopy
import cartopy.crs as ccrs
import pandas as pd
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=['c','#ff7f00', '#4daf4a',
                  '#f781bf', '#984ea3',
                                     '#e41a1c', '#dede00'])#r", "k", "c","m","darkorange",""]) 

def main():
    
    gaw  = pd.read_csv('lat_lons_GAW_VOCs.csv')
    #noaa = pd.read_csv('lat_lons_NOAA_VOCs.csv')
    na = pd.read_csv('lat-lon_n-america.csv')
    pa = pd.read_csv('lat-lon_pacific.csv')
    at = pd.read_csv('lat-lon_atlantic.csv')
    eu = pd.read_csv('lat-lon_europe.csv')
    As = pd.read_csv('lat-lon_asia.csv')
    sh = pd.read_csv('lat-lon_SH.csv')


    fig=plt.figure(figsize=(7,5))
    ax = fig.add_subplot(1,1,1, projection=ccrs.EqualEarth(), aspect='auto')
    ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', alpha=.2,facecolor='forestgreen'))
    ax.add_feature(cartopy.feature.NaturalEarthFeature('physical', 'ocean', '110m', edgecolor='face', alpha=.05,facecolor='navy'))

    ax.scatter( gaw.Lon, gaw.Lat, marker='^',s=100,edgecolor='k',c='gold',label='GAW site',transform=ccrs.PlateCarree(), zorder=200 )
    ax.scatter( na.lon, na.lat, marker='*',s=200,edgecolor='k',label='North America',transform=ccrs.PlateCarree(), zorder=100 )
    ax.scatter( eu.lon, eu.lat, marker='*',s=200,edgecolor='k',label='Europe/N. Africa',transform=ccrs.PlateCarree(), zorder=100 )
    ax.scatter( As.lon, As.lat, marker='*',s=200,edgecolor='k',label='Asia',transform=ccrs.PlateCarree(), zorder=100 )
    ax.scatter( at.lon, at.lat, marker='*',s=200,edgecolor='k',label='Atlantic',transform=ccrs.PlateCarree(), zorder=100 )
    ax.scatter( pa.lon, pa.lat, marker='*',s=200,edgecolor='k',label='Pacific',transform=ccrs.PlateCarree(), zorder=100 )
    ax.scatter( sh.lon, sh.lat, marker='*',s=200,edgecolor='k',label='Southern Hemisphere',  transform=ccrs.PlateCarree(), zorder=100 )

    ax.coastlines()
    ax.add_feature(cartopy.feature.BORDERS, alpha=.1)
    ax.set_extent([-179, 179.9, -89, 89])

    plt.legend(ncol=7,frameon=False, facecolor=None, loc=6, bbox_to_anchor=(0.032, -.05), fontsize=5.75)
    
    plt.tight_layout()
    plt.savefig( f'plots/site_map.png' , dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
