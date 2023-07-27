#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=diff_timeseries
#SBATCH --ntasks=1
#SBATCH --mem=48Gb
#SBATCH --partition=nodes
#SBATCH --time=00:45:00
#SBATCH --output=Logs/timeseries_%A.log
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
matplotlib.use('agg')
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
from CVAO_dict import CVAO_dict as d
import CVAO_tools as CV
plt.style.use('seaborn-darkgrid')

inputs=GC.get_arguments()
variable=inputs.var
if variable==None:
    variable='ethane'

k='2017regionsboth'
#rundirs=['fullchem_hourly_default', 'fullchem_2xEU_'+k, 'fullchem_2xNAM_'+k,'fullchem_2xAsia_'+k]
rundirs=['fullchem_hourly_default', '2xethane_EU', '2xethane_US','2xethane_Asia', '2xethane_Afr',
        '1.1xNOx_EU', '1.1xNOx_US','1.1xNOx_Asia', '1.1xNOx_Afr']#, '1.1xNOx_notSHP_Afr']
#rundirs=['fullchem_hourly_default', 'no_shipping_global','no_shipping_local']
string=''+k+'_sensitivity'

dfs=[]
for a in rundirs:
    print(a, '!!')
    var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='GEOS-Chem', year='2017')

    a_time=[]
    y = rp.find_nearest(lat, 16.9)
    x = rp.find_nearest(lon, -24.9)
    for t in range(len(atime)):
        v=var[t,0,y,x]
        a_time.append(v)
        atime[t]=atime[t].replace(hour=0, minute=0)

    a=pd.DataFrame({a:a_time}, index=atime)
    a=a.resample('D').mean()
    dfs.append(a)


diffs=[]
for n in range(1, len(dfs)):
    col=dfs[n].columns[0]
    diff = ( dfs[n].values - dfs[0].values )* 10.
    diff = diff[:,0]

    a=pd.DataFrame({col:diff}, index=dfs[0].index)
    diffs.append(a)

cv_df=pd.concat(diffs, axis=1, ignore_index=False)
print(cv_df)

labels=['2x EU C2H6','2x US C2H6','2x Asia C2H6','2x Africa C2H6','1.1x EU NOx','1.1x US NOx','1.1x Asia NOx', '1.1x Africa NOx']
#labels=rundirs[1:]#['2x EU C2H6','2x US C2H6','2x Asia C2H6','2x Africa C2H6','1.1x EU NOx','1.1x US NOx','1.1x Asia NOx', '1.1x Africa NOx']

colors=[ 'lightcoral','coral','red','darkred']
colors=['slategrey','royalblue','blue','midnightblue', 'lightcoral','coral','red','darkred']
styles=['-','dotted','dashed','dashdot','-','dotted','dashed','dashdot']
rp.timeseries_plot(diffs, variable, cvao_data=False, resample=False, string=string, 
            Ylabel='Cape Verde $O_3$ sensitivity (10*ppb)',labels=labels,colors=colors,styles=styles)

diffs=[]
for n in range(1, len(dfs)):
    col=dfs[n].columns[0]
    diff = - ( 100 - ( ( dfs[n].values / dfs[0].values ) * 100 ))
    diff = diff[:,0]
    a=pd.DataFrame({col:diff}, index=dfs[0].index)
    diffs.append(a)

rp.timeseries_plot(diffs, variable, cvao_data=False, resample=False, string='pc'+string, 
            Ylabel='CVAO $O_3$ sensitivity (%)',labels=labels,colors=colors,styles=styles)
