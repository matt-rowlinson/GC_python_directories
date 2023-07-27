#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=flasks_noaa
#SBATCH --ntasks=1
#SBATCH --mem=4Gb
#SBATCH --partition=interactive
#SBATCH --time=00:03:00
#SBATCH --output=Logs/flasks_noaa_%A.log
#SBATCH --open-mode=appendltruncate
import pandas as pd
import numpy as np
import sys
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib
plt.style.use('seaborn-darkgrid')
matplotlib.rc('font', family='P052') 
matplotlib.rc('font', serif='P052') 
matplotlib.rc('text', usetex='false') 
#matplotlib.rcParams.update({'font.size': 22})
from scipy import odr as odr

def func(p, x):
    a, b, = p
    return a * x + b

def calc_stats(obs, model):
    absError= model - obs
    SE = np.square(absError)
    MSE= np.mean(SE)
    RMSE=np.round( np.sqrt(MSE), 3)
    R2 = np.round( 1. - (np.var(absError) / np.var(obs)), 3)
    return RMSE, R2

def match_axis_lims(ax):
    a,b=ax.get_xlim()
    c,d=ax.get_ylim()
    ax.set_xlim(right=np.max([b,d]), left=np.min([a,c]))
    ax.set_ylim( ax.get_xlim() )
    return

def draw_1to1_line(ax):
    line = mlines.Line2D([0, 1], [0, 1], color='darkgrey', alpha=1., zorder=0)
    transform = ax.transAxes
    line.set_transform(transform)
    ax.add_line(line)
    return

def get_odr(x,y):
    from scipy import odr
    quad_model = odr.Model(func)

    data = odr.RealData(np.log10(np.squeeze(x)), np.log10(np.squeeze(y)))
    odr = odr.ODR(data, quad_model, beta0=[0., 1.])
    out = odr.run()
    popt = out.beta
    perr = out.sd_beta
    for i in range(len(popt)):
        nstd = 10. # to draw 5-sigma intervals
        popt_up = popt + nstd * perr
        popt_dw = popt - nstd * perr
        x_fit = np.linspace(min(np.log10(x)), max(np.log10(x)), len(x))
        fit = func(popt, x_fit)
    return x_fit, fit


def main():
    version='13.1.2'
    ## Read all the global data first
    ethane_obs    = pd.read_csv(f'v{version}_csvs/noaa_ethane_observations_2016.csv').values
    propane_obs   = pd.read_csv(f'v{version}_csvs/noaa_propane_observations_2016.csv').values
    A_alk4_obs= pd.read_csv(f'v{version}_csvs/ALK4_obs_{version}_DLim.csv', delimiter=',')#.values 
    alk4_obs  = A_alk4_obs[(A_alk4_obs > 0).all(1)].values
    ## Read all the NH data 
    NHethane_obs    = pd.read_csv(f'v{version}_csvs/NHnoaa_ethane_observations_2016.csv').values
    NHpropane_obs   = pd.read_csv(f'v{version}_csvs/NHnoaa_propane_observations_2016.csv').values
    B_alk4_obs  = pd.read_csv(f'v{version}_csvs/NHALK4_obs_{version}_DLim.csv', delimiter=',')#.values 
    NHalk4_obs  = B_alk4_obs[(B_alk4_obs > 0).all(1)].values
    ## Read all the SH data 
    SHethane_obs    = pd.read_csv(f'v{version}_csvs/SHnoaa_ethane_observations_2016.csv').values
    SHpropane_obs   = pd.read_csv(f'v{version}_csvs/SHnoaa_propane_observations_2016.csv').values
    C_alk4_obs  = pd.read_csv(f'v{version}_csvs/SHALK4_obs_{version}_DLim.csv', delimiter=',')#.values 
    SHalk4_obs  = C_alk4_obs[(C_alk4_obs > 0).all(1)].values
   
    version='13.1.2'
    ethane_ceds     = pd.read_csv(f'v{version}_csvs/ethane_ceds_only_{version}_2016.csv').values
    NHethane_ceds   = pd.read_csv(f'v{version}_csvs/NHethane_ceds_only_{version}_2016.csv').values
    SHethane_ceds   = pd.read_csv(f'v{version}_csvs/SHethane_ceds_only_{version}_2016.csv').values
    propane_ceds    = pd.read_csv(f'v{version}_csvs/propane_ceds_only_{version}_2016.csv').values
    NHpropane_ceds  = pd.read_csv(f'v{version}_csvs/NHpropane_ceds_only_{version}_2016.csv').values
    SHpropane_ceds  = pd.read_csv(f'v{version}_csvs/SHpropane_ceds_only_{version}_2016.csv').values
    alk4_ceds       = pd.read_csv(f'v{version}_csvs/ALK4_ceds_only_{version}_DLim.csv', delimiter=',')#.values *1e12
    alk4_ceds       = alk4_ceds[(A_alk4_obs > 0).all(1)].values 
    NHalk4_ceds     = pd.read_csv(f'v{version}_csvs/NHALK4_ceds_only_{version}_DLim.csv', delimiter=',')#.values *1e12
    NHalk4_ceds     = NHalk4_ceds[(B_alk4_obs > 0).all(1)].values 
    SHalk4_ceds     = pd.read_csv(f'v{version}_csvs/SHALK4_ceds_only_{version}_DLim.csv', delimiter=',')#.values *1e12
    SHalk4_ceds     = SHalk4_ceds[(C_alk4_obs > 0).all(1)].values

    version='13.4.0'
    ethane_geo     = pd.read_csv(f'v{version}_csvs/ethane_geo_only_{version}_2016.csv').values
    NHethane_geo   = pd.read_csv(f'v{version}_csvs/NHethane_geo_only_{version}_2016.csv').values
    SHethane_geo   = pd.read_csv(f'v{version}_csvs/SHethane_geo_only_{version}_2016.csv').values
    propane_geo    = pd.read_csv(f'v{version}_csvs/propane_geo_only_{version}_2016.csv').values
    NHpropane_geo  = pd.read_csv(f'v{version}_csvs/NHpropane_geo_only_{version}_2016.csv').values
    SHpropane_geo  = pd.read_csv(f'v{version}_csvs/SHpropane_geo_only_{version}_2016.csv').values
    alk4_geo       = pd.read_csv(f'v{version}_csvs/ALK4_geo_only_{version}_DLim.csv', delimiter=',')#.values *1e12
    alk4_geo       = alk4_geo[(A_alk4_obs > 0).all(1)].values 
    NHalk4_geo     = pd.read_csv(f'v{version}_csvs/NHALK4_geo_only_{version}_DLim.csv', delimiter=',')#.values *1e12
    NHalk4_geo     = NHalk4_geo[(B_alk4_obs > 0).all(1)].values 
    SHalk4_geo     = pd.read_csv(f'v{version}_csvs/SHALK4_geo_only_{version}_DLim.csv', delimiter=',')#.values *1e12
    SHalk4_geo     = SHalk4_geo[(C_alk4_obs > 0).all(1)].values 

    y1_plots= [ ethane_ceds, propane_ceds, alk4_ceds, 
                NHethane_ceds, NHpropane_ceds, NHalk4_ceds,
                SHethane_ceds, SHpropane_ceds, SHalk4_ceds ]
    y3_plots= [ ethane_geo, propane_geo, alk4_geo, 
                NHethane_geo, NHpropane_geo, NHalk4_geo,
                SHethane_geo, SHpropane_geo, SHalk4_geo ]
    x_plots= [ ethane_obs, propane_obs, alk4_obs,
               NHethane_obs, NHpropane_obs, NHalk4_obs,
               SHethane_obs, SHpropane_obs, SHalk4_obs ]

    ## Plot all data and odr
    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(12,12))
    axes=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
    cs = ['#0000ff','#ff0000', 'purple','y'] ; fn='P052'
    cs = ['b','y','#e41a1c']
    for n in range(len(x_plots)):
            print( n )

            Pdev  = axes[n].scatter(    x_plots[n], y1_plots[n], c=cs[0], alpha=.2, zorder=1, label='v13.1.2')
            x_fit1, fit1 = get_odr( x_plots[n], y1_plots[n] )
            axes[n].plot(10**x_fit1, 10**fit1, c='darkblue', lw=2, label='__best fit curve')
            
            
            Pdev3 = axes[n].scatter(    x_plots[n], y3_plots[n], c=cs[2], marker='o', alpha=.5, zorder=3, 
                                        label='v13.4.0' )
            x_fit1, fit1 = get_odr( x_plots[n], y1_plots[n] )
            x_fit3, fit3 = get_odr( x_plots[n], y3_plots[n] )

            axes[n].plot(10**x_fit1, 10**fit1, c='darkblue', lw=2, label='__best fit curve')
            axes[n].plot(10**x_fit3, 10**fit3, c='darkred',   lw=2, label='__best fit curve')
            #axes[n].set_yscale('log') ; axes[n].set_xscale('log')
            match_axis_lims(axes[n]) ; draw_1to1_line(axes[n])
    
    ax1.legend(fontsize=13)
    ax1.set_ylabel(f'Model ethane (ppt)')
    ax2.set_ylabel(f'Model propane (ppt)')
    ax3.set_ylabel(f'Model ALK4 (ppt)')
    ax4.set_ylabel(f'Model ethane (ppt)')
    ax5.set_ylabel(f'Model propane (ppt)')
    ax6.set_ylabel(f'Model ALK4 (ppt)')
    ax7.set_ylabel(f'Model ethane (ppt)')
    ax8.set_ylabel(f'Model propane (ppt)')
    ax9.set_ylabel(f'Model ALK4 (ppt)')
    ax7.set_xlabel(f'Observed ethane (ppt)')
    ax8.set_xlabel(f'Observed propane (ppt)')
    ax9.set_xlabel(f'Observed ALK4 (ppt)')
     
    ax2.set_title('Global', fontname=fn, fontsize=18, weight='bold')
    ax5.set_title('Northern Hemisphere', fontname=fn, fontsize=18, weight='bold')
    ax8.set_title('Southern Hemisphere', fontname=fn, fontsize=18, weight='bold')
    plt.tight_layout()
    plt.savefig(f'plots/bases-scale_figure.png')
    plt.close()


if __name__=="__main__":
    main()
