#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=c2h6_noaa
#SBATCH --ntasks=1
#SBATCH --mem=134Gb
#SBATCH --partition=nodes
#SBATCH --time=01:45:00
#SBATCH --output=Logs/c2h6_noaa_%A.log
#SBATCH --open-mode=appendltruncate
import pandas as pd
import glob
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap
import re
import random
import matplotlib.pyplot as plt
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
plt.style.use('seaborn-darkgrid')

import matplotlib.patches as mpatches
from matplotlib.patches import Polygon

def draw_screen_poly( lons, lats, m, c='red'):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( list(xy), edgecolor='k', facecolor=c, alpha=0.3)
    plt.gca().add_patch(poly)

def plot_sites_on_map(mapped):
    f,ax = plt.subplots(figsize=(12,6))
    m = Basemap()
    m.bluemarble()
    for i in range(len(mapped)):
        m.scatter(mapped[i][1],mapped[i][2],marker='*',color='y',s=150,zorder=10, label=mapped[i][0])
        plt.text(mapped[i][1]+random.uniform(-2,2), mapped[i][2]+random.uniform(-.5,.5),
                 mapped[i][0], weight='bold', c='y')
    plt.savefig('plots/noaa_sitemap.png')
    plt.close()
    
def main():

    m = rp.get_basemap( resolution='l',lsmask=True, lines=False, lllon=-90, urlon=35, lllat=-10, urlat=70)
    m.bluemarble()
    # CVAO 
    x = -24.4 ; y = 16.9
    m.scatter( x, y, marker='*', s=150, zorder=10, c='gold')

    CV_big_lats=[y-5, y+22, y+22, y-5]
    CV_big_lons=[x-10, x-10, x+7, x+7]
    draw_screen_poly(CV_big_lons, CV_big_lats, m, c='y')
    CV_broad= mpatches.Patch(color='y', alpha=.3, label='CV_broad')

    CV_N_lats=[y, y+8, y+8, y]
    CV_N_lons=[x-7, x-7, x+7, x+7]
    draw_screen_poly(CV_N_lons, CV_N_lats, m)
    CV_N= mpatches.Patch(color='red', alpha=.3, label='CV_North')

    CV_NE_lats=[y, y+20, y+20, y]
    CV_NE_lons=[x, x, x+7, x+7]
    draw_screen_poly(CV_NE_lons, CV_NE_lats, m, c='g')
    CV_NE= mpatches.Patch(color='g', alpha=.3, label='CV_NE')

    # Mace Head
    x = -9.9 ; y = 53.3
    m.scatter( x, y, marker='*', s=150, zorder=10, c='b')

    MH_bigE_lats=[y-10, y+10, y+10, y-10]
    MH_bigE_lons=[x-25, x-25, x, x]
    draw_screen_poly(MH_bigE_lons, MH_bigE_lats, m, c='skyblue')
    MH_East= mpatches.Patch(color='skyblue', alpha=.3, label='MH_East')

    MH_E_lats=[y-5, y+5, y+5, y-5]
    MH_E_lons=[x-15, x-15, x, x]
    draw_screen_poly(MH_E_lons, MH_E_lats, m, c='darkblue')
    MH_East_small= mpatches.Patch(color='darkblue', alpha=.3, label='MH_E_small')

    Atl_lats=[30, 60,60,30]
    Atl_lons=[-52, -52, -10, -10]
    draw_screen_poly(Atl_lons, Atl_lats, m, c='purple')
    Atl= mpatches.Patch(color='purple', alpha=.3, label='North Atlantic')

    plt.legend(handles=[CV_broad, CV_N, CV_NE, MH_East, MH_East_small, Atl], 
                    framealpha=1., facecolor='w',bbox_to_anchor=[.03,.2],frameon=True, loc='center left')
    plt.tight_layout()
    plt.savefig('plots/areas_map.png', dpi=400)
    plt.close()
    
if __name__=="__main__":
    main()
