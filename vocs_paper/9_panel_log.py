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
    ## Read all the global data first
    ethane_obs    = pd.read_csv('VOC_csvs/noaa_ethane_observations_2016.csv').values
    propane_obs   = pd.read_csv('VOC_csvs/noaa_propane_observations_2016.csv').values
    ethane_ceds   = pd.read_csv('VOC_csvs/ethane_ceds_only_13.4.0_2016.csv').values
    propane_ceds  = pd.read_csv('VOC_csvs/propane_ceds_only_13.4.0_2016.csv').values
    ethane_scale  = pd.read_csv('VOC_csvs/ethane_all_scale_13.4.0_2016.csv').values
    propane_scale = pd.read_csv('VOC_csvs/propane_all_scale_13.4.0_2016.csv').values   
    A_alk4_obs = pd.read_csv('csvs/ALK4_obs_v13.4.DLim.csv', delimiter=',')
    alk4_obs = A_alk4_obs[(A_alk4_obs > 0).all(1)].values
    alk4_ceds = pd.read_csv('csvs/ALK4_model2_v13.4.DLim.csv', delimiter=',')
    alk4_ceds = alk4_ceds[(A_alk4_obs > 0).all(1)].values #* 1e12
    alk4_scale= pd.read_csv('csvs/ALK4_MEIC_v13.4.DLim.csv', delimiter=',')
    alk4_scale = alk4_scale[(A_alk4_obs > 0).all(1)].values #* 1e12
   
    ## Read all the NH data 
    NHethane_obs    = pd.read_csv('VOC_csvs/NHnoaa_ethane_observations_2016.csv').values
    NHpropane_obs   = pd.read_csv('VOC_csvs/NHnoaa_propane_observations_2016.csv').values
    NHethane_ceds   = pd.read_csv('VOC_csvs/NHethane_ceds_only_13.4.0_2016.csv').values
    NHpropane_ceds  = pd.read_csv('VOC_csvs/NHpropane_ceds_only_13.4.0_2016.csv').values
    NHethane_scale  = pd.read_csv('VOC_csvs/NHethane_all_scale_13.4.0_2016.csv').values
    NHpropane_scale = pd.read_csv('VOC_csvs/NHpropane_all_scale_13.4.0_2016.csv').values   
    A_alk4_obs = pd.read_csv('csvs/NHALK4_obs_v13.4.DLim.csv', delimiter=',')#.values 
    NHalk4_obs = A_alk4_obs[(A_alk4_obs > 0).all(1)].values
    NHalk4_ceds = pd.read_csv('csvs/NHALK4_model2_v13.4.DLim.csv', delimiter=',')#.values *1e12
    NHalk4_ceds = NHalk4_ceds[(A_alk4_obs > 0).all(1)].values #* 1e12
    NHalk4_scale= pd.read_csv('csvs/NHALK4_MEIC_v13.4.DLim.csv', delimiter=',')#.values *1e12
    NHalk4_scale = NHalk4_scale[(A_alk4_obs > 0).all(1)].values #* 1e12

    ## Read all the SH data 
    SHethane_obs    = pd.read_csv('VOC_csvs/SHnoaa_ethane_observations_2016.csv').values
    SHpropane_obs   = pd.read_csv('VOC_csvs/SHnoaa_propane_observations_2016.csv').values
    SHethane_ceds   = pd.read_csv('VOC_csvs/SHethane_ceds_only_13.4.0_2016.csv').values
    SHpropane_ceds  = pd.read_csv('VOC_csvs/SHpropane_ceds_only_13.4.0_2016.csv').values
    SHethane_scale  = pd.read_csv('VOC_csvs/SHethane_all_scale_13.4.0_2016.csv').values
    SHpropane_scale = pd.read_csv('VOC_csvs/SHpropane_all_scale_13.4.0_2016.csv').values   
    A_alk4_obs = pd.read_csv('csvs/SHALK4_obs_v13.4.DLim.csv', delimiter=',')#.values 
    SHalk4_obs = A_alk4_obs[(A_alk4_obs > 0).all(1)].values
    SHalk4_ceds = pd.read_csv('csvs/SHALK4_model2_v13.4.DLim.csv', delimiter=',')#.values *1e12
    SHalk4_ceds = SHalk4_ceds[(A_alk4_obs > 0).all(1)].values #* 1e12
    SHalk4_scale= pd.read_csv('csvs/SHALK4_MEIC_v13.4.DLim.csv', delimiter=',')#.values *1e12
    SHalk4_scale = SHalk4_scale[(A_alk4_obs > 0).all(1)].values #* 1e12

    y1_plots= [ ethane_ceds, propane_ceds, alk4_ceds,
                NHethane_ceds, NHpropane_ceds, NHalk4_ceds,
                SHethane_ceds, SHpropane_ceds, SHalk4_ceds ]
    y2_plots= [ ethane_scale, propane_scale, alk4_scale, 
                NHethane_scale, NHpropane_scale, NHalk4_scale,
                SHethane_scale, SHpropane_scale, SHalk4_scale ]
    x_plots= [ ethane_obs, propane_obs, alk4_obs,
               NHethane_obs, NHpropane_obs, NHalk4_obs,
               SHethane_obs, SHpropane_obs, SHalk4_obs ]

    ## Plot all data and odr
    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(12,12))
    axes=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
    cs = ['#0000ff','#ff0000'] ; fn='P052'
    for n in range(len(x_plots)):
            #print( n )
            Pdev  = axes[n].scatter(    x_plots[n], y1_plots[n], c=cs[0], alpha=.2, zorder=1)
            Pdev2 = axes[n].scatter(    x_plots[n], y2_plots[n], c=cs[1], marker='+', alpha=.5, zorder=2 )
            x_fit1, fit1 = get_odr( x_plots[n], y1_plots[n] )
            x_fit2, fit2 = get_odr( x_plots[n], y2_plots[n] )
            axes[n].plot(10**x_fit1, 10**fit1, c='darkblue', lw=2, label='best fit curve')
            axes[n].plot(10**x_fit2, 10**fit2, c='darkred', lw=2, label='best fit curve')
            
            axes[n].set_yscale('log') ; axes[n].set_xscale('log')
            match_axis_lims(axes[n]) ; draw_1to1_line(axes[n])

    ax1.set_ylabel(f'Model ethane (ppt)')
    ax4.set_ylabel(f'Model ethane (ppt)')
    ax7.set_ylabel(f'Model ethane (ppt)')
    ax7.set_xlabel(f'Observed ethane (ppt)')
    ax8.set_xlabel(f'Observed propane (ppt)')
    ax9.set_xlabel(f'Observed ALK4 (ppt)')
     
    ax2.set_title('Global', fontname=fn, fontsize=18, weight='bold')
    ax5.set_title('Northern Hemisphere', fontname=fn, fontsize=18, weight='bold')
    ax8.set_title('Southern Hemisphere', fontname=fn, fontsize=18, weight='bold')
    plt.tight_layout()
    plt.savefig('plots/9_panel_figure_13.4.png')
    plt.close()

    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(12,12))
    axes=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9]
    for n in range(len(x_plots)):
            print( n )
            Pdev  = axes[n].scatter(    x_plots[n], y1_plots[n], c=cs[0], alpha=.2, zorder=1)
            x_fit1, fit1 = get_odr( x_plots[n], y1_plots[n] )
            x_fit2, fit2 = get_odr( x_plots[n], y2_plots[n] )
            axes[n].plot(10**x_fit1, 10**fit1, c='darkblue', lw=2, label='best fit curve')
            
            axes[n].set_yscale('log') ; axes[n].set_xscale('log')
            match_axis_lims(axes[n]) ; draw_1to1_line(axes[n])

    ax1.set_ylabel(f'Model ethane (ppt)') 
    ax4.set_ylabel(f'Model ethane (ppt)')
    ax7.set_ylabel(f'Model ethane (ppt)')  
    ax7.set_xlabel(f'Observed ethane (ppt)')
    ax8.set_xlabel(f'Observed propane (ppt)')
    ax9.set_xlabel(f'Observed ALK4 (ppt)')
    
    ax2.set_title('Global', fontname=fn, fontsize=18, weight='bold')
    ax5.set_title('Northern Hemisphere', fontname=fn, fontsize=18, weight='bold')
    ax8.set_title('Southern Hemisphere', fontname=fn, fontsize=18, weight='bold')
    plt.tight_layout()
    plt.savefig('plots/9_panel_figure_ced_only_13.4.png')
    plt.close()

    '''
    ax1.set_ylabel(f'Model ethane (ppt)')  ; ax1.set_xlabel(f'Observed ethane (ppt)')
    ax2.set_ylabel(f'Model propane (ppt)') ; ax2.set_xlabel(f'Observed propane (ppt)')
    ax3.set_ylabel(f'Model ALK4 (ppt)')    ; ax3.set_xlabel(f'Observed ALK4 (ppt)')
    ax4.set_ylabel(f'Model ethane (ppt)')  ; ax4.set_xlabel(f'Observed ethane (ppt)')
    ax5.set_ylabel(f'Model propane (ppt)') ; ax5.set_xlabel(f'Observed propane (ppt)')
    ax6.set_ylabel(f'Model ALK4 (ppt)')    ; ax6.set_xlabel(f'Observed ALK4 (ppt)')
    ax7.set_ylabel(f'Model ethane (ppt)')  ; ax7.set_xlabel(f'Observed ethane (ppt)')
    ax8.set_ylabel(f'Model propane (ppt)') ; ax8.set_xlabel(f'Observed propane (ppt)')
    ax9.set_ylabel(f'Model ALK4 (ppt)')    ; ax9.set_xlabel(f'Observed ALK4 (ppt)')
    '''

    return

if __name__=="__main__":
    main()
