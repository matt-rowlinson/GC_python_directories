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

from matplotlib.patches import Polygon
def draw_screen_poly( lons, lats, m):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( list(xy), edgecolor='k', facecolor='darkorange', alpha=0.3)
    plt.gca().add_patch(poly)

def plot_sites_on_map(mapped):
    f,ax = plt.subplots(figsize=(12,4))
    m = Basemap()
    #m.bluemarble()
    for i in range(len(mapped)):
        m.scatter(mapped[i][1],mapped[i][2],marker='*',color='y',s=150,zorder=10, label=mapped[i][0])
        plt.text(mapped[i][1]+random.uniform(-2,2), mapped[i][2]+random.uniform(-.5,.5),
                 mapped[i][0], weight='bold', c='y')
    plt.savefig('plots/noaa_sitemap.png')
    plt.close()
    
def main():
    nam_lats=[59.95, 20.05, 20.05, 59.95]
    nam_lons=[-139.95, -139.95, -50.05, -50.05]

    eur_lats=[30.05, 89.95, 89.95, 30.05]
    eur_lons=[-29.95, -29.95, 89.95, 89.95]

    asi_lats=[0.0, 50.0, 50.0, 0.0]
    asi_lons=[89.95, 89.95, 165.0, 165.0]
    
    f = plt.subplots( figsize=(12,6) )
    m = rp.get_basemap()
    #m = Basemap(projection='cyl',lon_0=0)
    #m.draw_coastlines()
    #m.bluemarble()

    draw_screen_poly(nam_lons, nam_lats, m)
    draw_screen_poly(eur_lons, eur_lats, m)
    draw_screen_poly(asi_lons, asi_lats, m)

    plt.tight_layout()
    plt.savefig('paper_figures/scaling_map.png', dpi=400)
    plt.close()
    
if __name__=="__main__":
    main()
