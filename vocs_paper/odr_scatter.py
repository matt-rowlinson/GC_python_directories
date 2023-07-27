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
plt.style.use('seaborn-darkgrid')
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
    ## Read all the data first
    ethane_obs    = pd.read_csv('csvs/noaa_ethane_observations_2016.csv').values
    propane_obs   = pd.read_csv('csvs/noaa_propane_observations_2016.csv').values

    ethane_ceds   = pd.read_csv('csvs/ethane_ceds_only_13.1.2_2016.csv').values
    propane_ceds  = pd.read_csv('csvs/propane_ceds_only_13.1.2_2016.csv').values

    ethane_scale  = pd.read_csv('csvs/ethane_asia-meic_scale_13.1.2_2016.csv').values
    propane_scale = pd.read_csv('csvs/propane_asia-meic_scale_13.1.2_2016.csv').values
    
    A_alk4_obs = pd.read_csv('csvs/ALK4_obs.csv', delimiter=',')#.values 
    alk4_obs = A_alk4_obs[(A_alk4_obs > 0).all(1)].values

    alk4_model2 = pd.read_csv('csvs/ALK4_model2.csv', delimiter=',')#.values *1e12
    alk4_model2 = alk4_model2[(A_alk4_obs > 0).all(1)].values * 1e12

    alk4_model3= pd.read_csv('csvs/ALK4_MEIC.csv', delimiter=',')#.values *1e12
    alk4_model3 = alk4_model3[(A_alk4_obs > 0).all(1)].values * 1e12

    ## Plot all data and odr
    fig, ((ax1,ax2,ax3)) = plt.subplots(1,3,figsize=(15,5))
    cs = ['#0000ff','#ff0000']
    cs = ['k','k']

    Pdev  = ax1.scatter( ethane_obs, ethane_ceds, c=cs[0], alpha=.2, zorder=1)
    Pdev2 = ax1.scatter( ethane_obs, ethane_scale, c=cs[1], marker='+', alpha=.5, zorder=2 )
    x_fit1, fit1 = get_odr( ethane_obs, ethane_ceds )
    x_fit2, fit2 = get_odr( ethane_obs, ethane_scale )
    ax1.plot(10**x_fit1, 10**fit1, c='darkblue', lw=2, label='best fit curve')
    ax1.plot(10**x_fit2, 10**fit2, c='darkred', lw=2, label='best fit curve')

    Pdev  = ax2.scatter( propane_obs, propane_ceds , c=cs[0], alpha=.2, zorder=1, label='Base')
    Pdev2 = ax2.scatter( propane_obs, propane_scale, c=cs[1], marker='+', alpha=.5, zorder=2, label='Scaled' )
    ax2.legend()
    x_fit1, fit1 = get_odr( propane_obs, propane_ceds )
    x_fit2, fit2 = get_odr( propane_obs, propane_scale )
    ax2.plot(10**x_fit1, 10**fit1, c='darkblue', lw=2, label='__best fit curve')
    ax2.plot(10**x_fit2, 10**fit2, c='darkred' , lw=2, label='__best fit curve')
    plt.legend()

    Pdev  = ax3.scatter( alk4_obs, alk4_model2 , c=cs[0], alpha=.05, zorder=1)
    Pdev2 = ax3.scatter( alk4_obs, alk4_model3 , c=cs[1], marker='+',  alpha=.15, zorder=2 )
    x_fit1, fit1 = get_odr( alk4_obs, alk4_model2 )
    x_fit2, fit2 = get_odr( alk4_obs, alk4_model3 )
    ax3.plot(10**x_fit1, 10**fit1, c='darkblue', lw=2, label='__best fit curve')
    ax3.plot(10**x_fit2, 10**fit2, c='darkred',  lw=2, label='__best fit curve')

    #ax1.set_yscale('log')
    #ax1.set_xscale('log')
    #ax2.set_yscale('log')
    #ax2.set_xscale('log')
    #ax3.set_yscale('log')
    #ax3.set_xscale('log')

    ax1.set_ylabel(f'Model ethane (ppt)')
    ax1.set_xlabel(f'Observed ethane (ppt)')
    ax2.set_ylabel(f'Model propane (ppt)')
    ax2.set_xlabel(f'Observed propane (ppt)')
    ax3.set_ylabel(f'Model ALK4 (ppt)')
    ax3.set_xlabel(f'Observed ALK4 (ppt)')
    
    match_axis_lims(ax1)
    match_axis_lims(ax2)
    match_axis_lims(ax3)

    draw_1to1_line(ax1)
    draw_1to1_line(ax2)
    draw_1to1_line(ax3)

    ax3.set_xlim( left=0., right=1850)
    ax3.set_ylim( ax3.get_xlim() ) 
    
    plt.tight_layout()
    plt.savefig('plots/TEST.figure_5.png')
    plt.close()

    return

if __name__=="__main__":
    main()
