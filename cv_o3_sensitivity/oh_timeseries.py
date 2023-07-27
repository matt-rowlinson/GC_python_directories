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

rundirs=['fullchem_hourly_default', 'fullchem_2xEU_ethane', 'fullchem_2xNAM_ethane','fullchem_2xAsia_ethane']
rundirs=['fullchem_hourly_default', 'fullchem_onlyCEDS','fullchem_onlyTzompa']
string='onlyCEDS_onlyTzompa_2016'

dfs=[]
for a in rundirs:
    var, lat, lon, lev, atime = GC.get_gc_var(rundir=a, variable=variable, collection='ConcAfterChem',
                                            version='GEOS-Chem', year='2016')
    
    a_time=[] 
    y = rp.find_nearest(lat, 16.9)
    x = rp.find_nearest(lon, -24.9)
    for t in range(len(atime)):
        v=var[t,0,y,x]
        a_time.append(v)
    print(atime)
    a=pd.DataFrame({a:a_time}, index=atime)
    dfs.append(a)

print(dfs)
cv_df=pd.concat(dfs, axis=1, ignore_index=False)
#print(cv_df)

dfs[0]=dfs[0].resample('D').mean()

labels=['GC v13.0.1', '2x EU ethane','2x US ethane','2x Asian ethane']
labels=['HEMCO_default','only_CEDS', 'only_Tzompa']
styles=['dotted','dashed','dashdot','dotted','dashed','dashdot']
rp.timeseries_plot(dfs, variable, cvao_data=False,  resample='D', string=string, labels=labels, styles=styles)
rp.timeseries_plot(dfs, variable, cvao_data=False, std=True, resample='D', string=string, labels=labels, styles=styles)
