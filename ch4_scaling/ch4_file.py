import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import sys
sys.path.append('/users/mjr583/python_lib/')
import RowPy as rp
import GC_tools as GC
import CVAO_tools as CV
from CVAO_dict import CVAO_dict as d


## Scaling file for GC
f='monthly.gridded.surface.methane.1979-2020.1x1.nc'
fh=Dataset(f)
lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
time=fh.variables['time']
ch4=fh.variables['SFC_CH4'][:]

## CH4 at Cape Verde
ladx=rp.find_nearest(lat, 16.9)
lodx=rp.find_nearest(lon, -24.9)
time=time[:]
time=rp.hours_to_datetime(time,start=1979)
cv_ch4=ch4[:,ladx,lodx]
cv_ch4 = pd.DataFrame({'CV':cv_ch4}, index=time)

## CH4 at Tenerife
ladx=rp.find_nearest(lat, 28.309)
lodx=rp.find_nearest(lon, -16.499)
ten_ch4=ch4[:,ladx,lodx]
ten_ch4 = pd.DataFrame({'Ten_CH4':ten_ch4}, index=time)

## CH4 at Tenerife
ladx=rp.find_nearest(lat, 18.9841)
lodx=rp.find_nearest(lon, -97.311)
mex_ch4=ch4[:,ladx,lodx]
mex_ch4 = pd.DataFrame({'Mex_CH4':mex_ch4}, index=time)


## CH4 from NOOA obspack
fh=Dataset('/users/mjr583/scratch/noaa_ch4/data/nc/ch4_izo_surface-flask_1_representative.nc')
ch4=fh.variables['value'][:]*1e9
time=fh.variables['time'][:]
new_time=rp.seconds_to_datetime(time,start=1970)
noaa=pd.DataFrame({'Value':ch4}, index=new_time)
noaa=noaa['2007':]
noaa_mon=noaa.resample('M').mean()

fh=Dataset('/users/mjr583/scratch/noaa_ch4/data/nc/ch4_mex_surface-flask_1_representative.nc')
ch4=fh.variables['value'][:]*1e9
time=fh.variables['time'][:]
new_time=rp.seconds_to_datetime(time,start=1970)
mex=pd.DataFrame({'Value':ch4}, index=new_time)
mex=mex['2007':]
mex_mon=mex.resample('M').mean()

## CH4 in 4x5 GEOS-Chem run at CVO and Tenerife and Mexico
rundir='global_4x5'
var, lat, lon, lev, time = GC.get_gc_var(rundir, 'CH4')
y = rp.find_nearest(lat, 16.9)
x = rp.find_nearest(lon, -24.9)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
gc_cv=pd.DataFrame({'Value':var_time}, index=time)

y = rp.find_nearest(lat, 28.209)
x = rp.find_nearest(lon, -16.499)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
gc_ten=pd.DataFrame({'Value':var_time}, index=time)

y = rp.find_nearest(lat, 18.9841)
x = rp.find_nearest(lon, -97.311)
var_time=[]
for t in range(len(time)):
    v=var[t,0,y,x]
    var_time.append(v)
gc_mex=pd.DataFrame({'Value':var_time}, index=time)


## CH4 at CVAO
cv=CV.get_from_merge(d['CH4'])
cv=cv['2007':]
cv=cv.resample('M').mean()


## Now just plot
plt.plot(cv_ch4['2007':], label='Scaling values at CVAO')
plt.plot(ten_ch4['2007':], label='Scaling values at Tenerife')
plt.plot(gc_cv, label='GC CVAO')
plt.plot(gc_ten, label='GC Tenerife')
plt.plot(cv, label='CVAO')
plt.plot(noaa_mon, label='Tenerife')
plt.legend()
plt.ylabel('$CH_4$ (ppbv)')
plt.savefig('./plots/ch4_all.png')
plt.close()

plt.plot(cv_ch4['2007':], label='Scaling values at CVAO')
plt.plot(ten_ch4['2007':], label='Scaling values at Tenerife')
plt.legend()
plt.ylabel('$CH_4$ (ppbv)')
plt.savefig('./plots/ch4_scalingonly.png')
plt.close()

plt.plot(noaa_mon, label='Surface obs at Tenerife', color='r')
plt.plot(mex_mon, label='Surface obs at Mexico', color='g')
plt.plot(cv, label='Surface obs at CVAO', color='b')
plt.legend()
plt.ylabel('$CH_4$ (ppbv)')
plt.savefig('./plots/ch4_obsonly.png')
plt.close()


plt.plot(ten_ch4['2007':], label='GC scaling values at Tenerife')
plt.plot(gc_ten, label='GEOS-Chem Tenerife')
plt.plot(noaa_mon, label='Surface obs at Tenerife', color='k')
plt.legend()
plt.ylabel('$CH_4$ (ppbv)')
plt.savefig('./plots/ch4_tenerife.png')
plt.close()

plt.plot(mex_ch4['2007':], label='GC scaling values at Mexico site')
plt.plot(gc_mex, label='GEOS-Chem Mexico site')
plt.plot(mex_mon, label='Surface obs at Mexico site', color='k')
plt.legend()
plt.ylabel('$CH_4$ (ppbv)')
plt.savefig('./plots/ch4_mexico.png')
plt.close()


plt.plot(cv_ch4['2007':], label='GC scaling values at CVAO')
plt.plot(gc_cv, label='GEOS-Chem CVAO')
plt.plot(cv, label='Surface obs at CVAO', color='k')
plt.legend()
plt.ylabel('$CH_4$ (ppbv)')
plt.savefig('./plots/ch4_capeverde.png')
plt.close()
