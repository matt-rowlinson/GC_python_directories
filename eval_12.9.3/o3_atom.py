import numpy as np
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime
import matplotlib.pyplot as plt
import glob
import sys
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import RowPy as rp
path='/users/mjr583/scratch/ATom/ATom_merge_1581/data/'

plt.style.use('seaborn-darkgrid')

altitude=[] ; timex=[] ; variable=[] ; var='Not yet'
for infile in sorted(glob.glob('%sMER-SAGA*.nc' %path)):
    print(infile)
    fh=Dataset(infile)
    try:
        data=fh.groups['UCATS-O3'].variables['O3_UO3'][:]
    except:
        print('O3 not present in %s' %infile)
        continue

    times = fh.variables['time'][:]
    t = datetime(2016,1,1,0,0) ## time is given in seconds since 2016-01-01. Converting into datetime object
    ts = (t-datetime(1970,1,1,0,0)).total_seconds()
    timex=[] 
    for i in range(len(times)):
        timex.append(datetime.fromtimestamp(times[i]+ts).strftime("%d/%m/%Y %H:%M:%S"))
    dates_list= [datetime.strptime(time, "%d/%m/%Y %H:%M:%S").strftime('"%d/%m/%Y %H:%M:%S"') for time in timex]

    var = pd.DataFrame({'date':dates_list, 'O3':data[:]})
    var = var.set_index(['date'])
    var.index = pd.to_datetime(var.index,format='"%d/%m/%Y %H:%M:%S"')

    alt = fh.groups['MMS'].variables['G_ALT'][:] * 1e-3
    var['altitude']=alt
    lon = fh.groups['MMS'].variables['G_LONG'][:]    
    var['longitude']=lon
    lat = fh.groups['MMS'].variables['G_LAT'][:]
    var['latitude']=lat
    
    variable.append(var)

#saga1=variable[0]
#saga2=variable[1]
#saga3=variable[2]
#df=pd.concat(variable)

x=-24.9
y=16.9

variable2=[]
for df in variable:
    df=df.loc[df['latitude'] > y-10]
    df=df.loc[df['latitude'] < y+10]
    df=df.loc[df['longitude'] > x-10]
    df=df.loc[df['longitude'] < x+10]
    variable2.append(df)

saga1=variable2[0]
saga2=variable2[1]
saga3=variable2[2]

print(saga3)

var, lat, lon, lev, time = GC.get_gc_var('tropchem_merra_4x5','O3', '12.9.3', year='2016')
lev=pd.read_csv('../info_files/GC_vertical_levels.csv')['Altitude (km)']
print(lev)
y = rp.find_nearest(lat, y)
x = rp.find_nearest(lon, x)

gc1=var[7,:,y,x]
gc2=var[1,:,y,x]
gc3=var[9,:,y,x]

plt.scatter(saga1.O3, saga1.altitude, label='ATom 1 (Aug)')#,color'r')
plt.scatter(saga2.O3, saga2.altitude, label='ATom 2 (Feb)')#, 'b')
plt.scatter(saga3.O3, saga3.altitude, label='ATom 3 (Oct)')#, 'g')
plt.plot(gc1[:34], lev[:34], label='Geos-Chem')
plt.plot(gc2[:34], lev[:34])
plt.plot(gc3[:34], lev[:34])

plt.xlabel('$O_3$ (ppbv)')
plt.ylabel('Altitude (km)')
plt.legend()
plt.xlim([10,100])
plt.savefig('plots/vertical_o3_withATom.png')
plt.close()
sys.exit()
