#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import glob
import numpy as np
import re
from netCDF4 import Dataset
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV

from compare_NMVOCs import v, v20, v18
v_voc   = np.array(v)
v20_voc = np.array(v20)
v18_voc = np.array(v18)

from compare_CEDS_versions import v2, v20, v18, x
v_ethane   = np.array(v2)
v20_ethane = np.array(v20)
v18_ethane = np.array(v18)
print(x)

print('Funky plot')
fig,ax = plt.subplots(figsize=(9,5))
fig.subplots_adjust(right=0.75)

ax1 = ax.twinx()
ax2 = ax.twinx()
ax2.spines['right'].set_position(("axes", 1.12))


## Dummy plot just for the legends
L1a, = ax.plot(  x[:len(v18_voc)], v18_voc ,c='k', ls='-', label='CEDS v2018 (to 2014)' )
L1b, = ax.plot(  x[:len(v20_voc)], v20_voc ,c='k',ls='--', label='CEDS v2020 (to 2017)' )
L1c, = ax.plot(  x, v_voc , ls=':',c='k', label='CEDS v2021 (to 2019)' )
ax.legend(handles=[L1a,L1b,L1c])

L2a, = ax.plot(  x[:len(v18_voc)], v18_voc ,c='b', ls='-', label='Total NMVOCs')
L2b, = ax1.plot(  x[:len(v18_ethane)], v18_ethane ,c='r',ls='-', label='Ethane')
L2c, = ax2.plot(  x[:len(v18_ethane)], v18_ethane / v18_voc , ls='-',c='g', label='Ethane / NMVOCs ratio')
ax1.legend(handles=[L2a,L2b,L2c])


## Real plots
p1a, = ax.plot(  x[:len(v18_voc)], v18_voc ,c='b', ls='-', label='CEDS v2018 (to 2014)' )
p1b, = ax.plot(  x[:len(v20_voc)], v20_voc ,c='b',ls='--', label='CEDS v2020 (to 2017)' )
p1c, = ax.plot(  x, v_voc , ls=':',c='b', label='CEDS v2021 (to 2019)' )

p2a, = ax1.plot(  x[:len(v18_ethane)], v18_ethane ,c='r',ls='-', label='CEDS v2018 (to 2014)' )
p2b, = ax1.plot(  x[:len(v20_ethane)], v20_ethane ,c='r',ls='--', label='CEDS v2020 (to 2017)' )
p2c, = ax1.plot(  x, v_ethane ,c='r',ls=':',  label='CEDS v2021 (to 2019)' )

p3a, = ax2.plot(  x[:len(v18_ethane)], v18_ethane / v18_voc, c='g', ls='-', label='CEDS v2018 (to 2014)' )
p3b, = ax2.plot(  x[:len(v20_ethane)], v20_ethane / v20_voc, c='g', ls='--', label='CEDS v2020 (to 2017)' )
p3c, = ax2.plot(  x, v_ethane / v_voc, ls=':', c='g', label='CEDS v2021 (to 2019)' )

ax.yaxis.label.set_color(p1a.get_color())
ax1.yaxis.label.set_color(p2a.get_color())
ax2.yaxis.label.set_color(p3a.get_color())

ax.set_ylabel( 'NMVOCs (Tg $yr^{-1}$)'  )
ax1.set_ylabel( 'C2H6 (Tg $yr^{-1}$)'   )
ax2.set_ylabel( 'C2H6 / NMVOCs Ratio'  )

tkw = dict(size=4, width=1.5)
ax.tick_params(axis='y', colors=p1a.get_color(), **tkw)
ax1.tick_params(axis='y', colors=p2a.get_color(), **tkw)
ax2.tick_params(axis='y', colors=p3a.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)

plt.savefig( 'plots/CEDS_versions_NMVOCs_ethane_ratio.png' )
plt.close()

sys.exit()

