#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=NOxtimeseries
#SBATCH --ntasks=1
#SBATCH --mem=18Gb
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

rundirs=['fullchem_hourly_default', '2xethane_EU', '2xethane_US','2xethane_Asia', '2xethane_Afr']
rundirs=['fullchem_hourly_default', 'no_shipping_global','no_shipping_local']

string='no_shipping_'

dfs=[]
for a in rundirs:
    print(a)
    var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, version='GEOS-Chem', year='2017')
    
    a_time=[] 
    y = rp.find_nearest(lat, 16.9)
    x = rp.find_nearest(lon, -24.9)
    for t in range(len(atime)):
        v=var[t,0,y,x]
        a_time.append(v)
    a=pd.DataFrame({a:a_time}, index=atime)
    dfs.append(a)

#print(dfs)
cv_df=pd.concat(dfs, axis=1, ignore_index=False)
#print(cv_df)

dfs[0]=dfs[0].resample('D').mean()

#labels=['GC v13.0.1', '1.5x Global Shipping','No shipping global']#2x US ethane','2x Asian ethane']
labels=['HEMCO_default','CEDS + Geological ethane', 'CEDS + Geological ethane + NEI']
labels=['HEMCO default', 'No Global Shipping', 'No Local Shipping']
styles=['dotted','dashed','dashdot','dotted','dashed','dashdot']
rp.timeseries_plot(dfs, variable, cvao_data=True,  string=string, labels=labels, styles=styles, logscale=True)
rp.timeseries_plot(dfs, variable, cvao_data=True, std=True, resample='D', string=string, labels=labels, styles=styles, logscale=True)
