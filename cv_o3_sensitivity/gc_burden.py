import sys
import pandas as pd
import numpy as np
import glob
from netCDF4 import Dataset
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
rundir='fullchem_hourly_default'

zboltz=1.3807e-23       #Boltzmann constant
Rd=287.05               #gas constant
avc=6.022e23            #Avogadro's constant
mm_da=avc*zboltz/Rd     #molar mass of dry air (approx 28 g/mol)
mm_o3=0.0479982
mm_moh=0.03204
mm_co=0.02801
mm_ch4 = 0.01604

var, lat,lon,lev,haltime = GC.get_gc_var(rundir='fullchem_noc2h6_ems', variable=variable, version='GEOS-Chem',year='2015')
VAR= np.mean(var,0) / float(d[variable]['scale'])

def get_tropopause_ps(rundir, version, year):
    for i,infile in enumerate(sorted(glob.glob('/users/mjr583/scratch/GC/%s/rundirs/%s/OutputDir/GEOSChem.StateMet.*%s*.nc4' % (version, rundir, year)))):
        print(infile)
        fh=Dataset(infile)
        Tlev=fh.variables['Met_TropLev'][:]
        ad=fh.variables['Met_AD'][:]
        return Tlev, ad

zboltz=1.3807e-23       #Boltzmann constantn
Rd=287.05               #gas constant
avc=6.022e23            #Avogadro's constant
mm_da=avc*zboltz/Rd     #molar mass of dry air (approx 28 g/mol)
mm_o3=0.0479982
mm_moh=0.03204
mm_co=0.02801
mm_ch4 = 0.01604
mm_c2h6= 0.03007

##---------------------MAIN SCRIPT---------------------------------------------##

## Get your Geos-Chem species variable. Need to be 3D (lev, lat, lon) so choose a time or average time. 
## Variable also needs to be in mol mol-1, NOT ppb/ppt

# This function gets the GC tropospause level (Tlev) and air mass per grid box (AD)
# You'll need to change the path in the function above
T, AD = get_tropopause_ps(rundir='fullchem_hourly_default', version='GEOS-Chem', year='2016')
T = T[0]
AD= AD[0]

for ilon in range(len(lon)):
    #print(ilon)
    for ilat in range(len(lat)):
        for ilev in range(1,48):            
            if (ilev > T[ilat,ilon]):
                VAR[ilev-1,ilat,ilon] = np.nan

mm=mm_c2h6
var=np.zeros(47) * np.nan
for i in range(47):
    var[i] = np.nansum( VAR[i,:,:] * AD[i,:,:] * mm/mm_da) * 1e-9

burd = np.round(np.nansum(var),2)
print('Tropospheric %s burden is %s Tg' %(variable, burd ) )
