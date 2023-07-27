#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=hemco_timeseries
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --partition=nodes
#SBATCH --time=00:00:10
#SBATCH --output=Logs/hemco_timeseries_%A.log
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import os
import numpy as np
import sys
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp

def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

days=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
sectors = ['Anthro','BioBurn']#,'Aircraft','Ship']

regions={
    False : False,
    'Asia' : { 'minlon' : 70, 'maxlon' : 165, 'minlat' : 0., 'maxlat' : 50. },
    'US'   : { 'minlon' : -160, 'maxlon' : -60, 'minlat' : 15., 'maxlat' : 72. },
    'EU'   : { 'minlon' : -20, 'maxlon' : 40, 'minlat' : 36., 'maxlat' : 72. },
    'NH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : 0., 'maxlat' : 90. },
    'SH'   : { 'minlon' : -180, 'maxlon' : 180, 'minlat' : -90., 'maxlat' : 0. },
        }
rundirs = ['base_run_1980','base_run_1990','base_run_2000','base_run_2005','base_run_2010','base_run_2015','base_run_2020']
species=['NO','CO','HNO2','C2H6','NO']
for speci in species:
    region = False
    #for rundir in rundirs:
    #    print(rundir)
    for region in regions:
        for rundir in rundirs:
            path='/users/mjr583/scratch/GC/13.1.2/rundirs/%s/OutputDir/' %rundir
            flist=find_file_list(path, ['HEMCO_diagnostics.'])
            a = xr.open_mfdataset(flist)

            if region:
                out=region
                minlon,maxlon,minlat,maxlat=(regions[region]['minlon'],regions[region]['maxlon'], 
                                             regions[region]['minlat'], regions[region]['maxlat'])
                a = a.\
                        where(a.lon>= minlon, drop=True).\
                        where(a.lon<= maxlon, drop=True).\
                        where(a.lat>= minlat, drop=True).\
                        where(a.lat<= maxlat, drop=True)
            else:
                out='global'
            
            AREA=a.AREA
            #AREA = rp.surface_area_earth(46,72)
            #print(AREA.values.sum())
            em=a.sum(dim='lev', keep_attrs=True)
            em = em[f'Emis{speci}_BioBurn']
            em=np.array(em)
            print(em.shape)
            days=np.array([31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31])
            new_em=[]
            
            for n in range(len(em)):
                x = em[n] * AREA[n] * ( 3600 * 24 * days[n]) * 1e-9
                new_em.append(x)
            em = np.array(new_em)
            print(rundir, out+': ', np.round(em.sum(),2), 'Tg ')

            a=a.sum(dim='lat', keep_attrs=True)
            a=a.sum(dim='lon', keep_attrs=True)
            a=a.sum(dim='lev', keep_attrs=True)
            
            a1d=False#np.zeros((12))
            #for i in sectors:
            #    #try:
            #    if a1d is False:
            try:
                a1d = a[f'Emis{speci}_BioBurn']
            except:
                continue
            #    else:
            #        try:
            #            if 'Geo' in i:
            #                print('here')
            #                a1d += (a[f'Emis{speci}_{i}' ] * 12)
            #                print('and here')
            #            else:
            #                 a1d += a[f'Emis{speci}_{i}' ] 
            #        except:
            #            continue


                #    print(a1d.sum())
                    #print(a1d.shape)
                #except:
                #    continue
            
               # if rundir=='fullchem_hourly_default':
               #     rundir='hemco_default'
               # else:
               #     label='only_ceds'
                    
                    #a1d.plot.step("-", label=rundir)
            plt.step(a['time'],a1d, label=rundir)
        plt.legend()
        plt.savefig(f'plots/onlyEMS_{speci}_{region}_BioBurn.png')
        plt.close()
        #sys.exit()
