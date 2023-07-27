#!/usr/bin/env python
"""
Six Panel Comparison Plots
--------------------------------------
This example script demonstrates the comparitive plotting capabilities of GCPy,
including single level plots as well as global zonal mean plots.
These comparison plots are frequently used to evaluate results from different runs / versions
of GEOS-Chem, but can also be used to compare results from different points in one run that
are stored in separate xarray datasets.
The example data described here is in lat/lon format, but the same code works equally
well for cubed-sphere (GCHP) data. 
"""

#xarray allows us to read in any NetCDF file, the format of most GEOS-Chem diagnostics,
#as an xarray Dataset
import xarray as xr
import numpy as np
import gcpy.plot as gcplot
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import math

#------------------------------------

refdir='/users/hc2134/scratch/GC/rundirs/merra2_4x5_standard_BryBudgets_Clybudgets_Iybudget_tagged'
devdir='/users/hc2134/scratch/GC/rundirs/merra2_4x5_HBrReversible_Brief_SSHOBrToBr2_SSHOBrToBrCl_SSHOClToCl2_SSHOClToBrCl_BrClToBr2_Cl2ToBrClandICl_Br2ToIBr_HOI_BryBudgets_ClyBudgets_IyBudget_tagged'

dev_ds = xr.open_dataset(devdir+'/OutputDir/GEOSChem.ProdLoss.20190701_0000z.nc4')
ref_ds = xr.open_dataset(refdir+'/OutputDir/GEOSChem.ProdLoss.20190701_0000z.nc4')


dev_concds = xr.open_dataset(devdir+'/OutputDir/GEOSChem.SpeciesConc.20190701_0000z.nc4')
ref_concds = xr.open_dataset(refdir+'/OutputDir/GEOSChem.SpeciesConc.20190701_0000z.nc4')


dev_metds = xr.open_dataset(devdir+'/OutputDir/GEOSChem.StateMet.20190701_0000z.nc4')
ref_metds = xr.open_dataset(refdir+'/OutputDir/GEOSChem.StateMet.20190701_0000z.nc4')

dev_speciesinfo=devdir+'/species_database.yml'
ref_speciesinfo=refdir+'/species_database.yml'


#--------
#exit()




#---------------------------

ref_ds['Prod_THOBrPFromBrClSS']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_THOBrPFromBrClSS.attrs["units"] = "molec cm-3 s-1"


ref_ds['Prod_THOBrPFromBr2SS']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_THOBrPFromBr2SS.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_THOClPFromCl2SS']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_THOClPFromCl2SS.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_THOIPFromIClSS']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_THOIPFromIClSS.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_THOIPFromIBrSS']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_THOIPFromIBrSS.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_THBrPFromSSBrm']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_THBrPFromSSBrm.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TIBrPFromICl']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TIBrPFromICl.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TI2PFromICl']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TI2PFromICl.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TI2PFromIBr']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TI2PFromIBr.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TIBrPFromBr2']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TIBrPFromBr2.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TBrClPFromCl2']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TBrClPFromCl2.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TBr2PFromBrCl']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TBr2PFromBrCl.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TIClPFromCl2SS']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TIClPFromCl2SS.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TI2PFromHOI']=ref_ds['Prod_TBrNO2PFromBr'] - ref_ds['Prod_TBrNO2PFromBr']
ref_ds.Prod_TI2PFromHOI.attrs["units"] = "molec cm-3 s-1"
ilev= 0 #30 #21

varlistBr = ['Prod_HBrToSS','Prod_TBrNO2PFromBr','Prod_THBrPFromSSBrm','Prod_THBrPFromBr',
	   'Prod_THBrPFromHOBr','Prod_TIBrPFromICl', 'Prod_TIBrPFromHOI',
	   'Prod_TIBrPFromIONO','Prod_TIBrPFromBr2', 'Prod_TIBrPFromIONO2',
           'Prod_THOBrPFromBrNO3', 'Prod_THOBrPFromO3', 'Prod_THOBrPFromBrClSS',
	   'Prod_THOBrPFromBr2wOH','Prod_THOBrPFromBrO', 'Prod_THOBrPFromBr2SS', 
           'Prod_TBrPFromHBr', 'Prod_TBrPFromBr2', 'Prod_TBrPFromBrNO3', 
           'Prod_TBrPFromHOBr', 'Prod_TBrPFromIBr', 'Prod_TBrPFromBrCl', 
	   'Prod_TBrPFromBrO', 'Prod_TBrPOTHPTW', 'Prod_TBrPFromBrNO2', 
	   'Prod_TBrClPFromClNO3', 'Prod_TBrClPFromCl2', 'Prod_TBrClPFromHOBr', 
           'Prod_TBrClPFromHOCl', 'Prod_TBrClPFromClNO2', 'Prod_TBrClPFromBrO', 
	   'Prod_TBrClPFromBrNO3', 'Prod_TBr2PFromBr', 'Prod_TBr2PFromHOBr', 
           'Prod_TBr2PFromBrCl', 'Prod_TBr2PFromBrO', 'Prod_TBrOPOTHPTW', 
	   'Prod_TBrOPFromHOBr', 'Prod_TBrOPFromHBr', 'Prod_TBrOPFromBr', 
	   'Prod_TBrOPFromBrNO3', 'Prod_TBrNO3PFromBrO', 'Prod_HOBr', 
           'Prod_Br', 'Prod_BrO', 'Prod_Br2', 
           'Prod_BrNO3', 'Prod_BrNO2', 'Prod_BrCl', 
	   'Prod_IBr', 'Prod_HBr'] 


varlistCl = ['Prod_TClPFromHOCl','Prod_TClPFromBrCl','Prod_TClPFromHCl',
                     'Prod_TClPFromClNO3','Prod_TClPFromClNO2','Prod_TClPFromCl2',
		     'Prod_TClPFromClOx','Prod_TClPOTHPTW','Prod_TClPFromICl',
		     'Prod_TCl2PFromOH','Prod_TCl2PFromClNO2','Prod_TCl2PFromClNO3',
		     'Prod_TCl2PFromHOCl','Prod_TCl2PFromClOx','Prod_THOClPFromClNO3',
		     'Prod_THOClPFromClOx','Prod_THOClPFromClNO2','Prod_THOClPFromCl2wOH',
		     'Prod_THOClPFromCl2SS',
		     'Prod_TClNO2PFromN2O5','Prod_THClPFromCl','Prod_THClPFromClOx',
		     'Prod_THClPFromHOCl','Prod_TIClPFromIONO','Prod_TIClPFromIONO2',
		     'Prod_TIClPFromHOI','Prod_TIClPFromClOx','Prod_TIClPFromCl2SS',
		     'Prod_TClNO3PFromClOx',
		     'Prod_TClOxPFromCl2','Prod_TClOxPOTHPTW','Prod_TClOxPFromHCl',
		     'Prod_TClOxPFromCl','Prod_TClOxPFromClNO3','Prod_TClOxPFromHOCl',
           'Prod_Cl', 'Prod_Cl2', 
           'Prod_ClNO3', 'Prod_ClNO2', 
	   'Prod_ICl', 'Prod_HCl'] 

varlistI = ['Prod_TIPFromHOI','Prod_TIPFromIBr','Prod_TIPFromICl',
                     'Prod_TIPFromHI','Prod_TIPFromINOx','Prod_TIPFromIONO2',
		     'Prod_TIPFromIOx','Prod_TIPFromI2','Prod_TIPOTHPTW',
		     'Prod_TINOxPFromI','Prod_TIONO2PFromIOx','Prod_TIONO2PFromI2',
		     'Prod_THOIPFromIOx','Prod_THOIPFromIClSS','Prod_THOIPFromIBrSS',
		     'Prod_THOIPFromI2','Prod_THIPFromI','Prod_TI2PFromICl',
		     'Prod_TI2PFrombiINOx','Prod_TI2PFromI','Prod_TI2PFromHOI',
		     'Prod_TI2PFromIBr','Prod_TIOxPFromHOI','Prod_TIOxPFromIONO2',
		     'Prod_TIOxPFromI','Prod_I','Prod_I2','Prod_IONO2','Prod_INOx',
		     'Prod_HI']






MW_g_Br =  np.zeros(len(varlistBr))
MW_g_Br[:]=79.90

MW_g_Cl =  np.zeros(len(varlistCl))
MW_g_Cl[:]=35.45

MW_g_I =  np.zeros(len(varlistI))
MW_g_I[:]=126.90


varlist=varlistBr# + varlistCl + varlistI 
#MW_g =np.concatenate(( MW_g_Br, MW_g_Cl, MW_g_I))
MW_g = MW_g_Br

#print(MW_g)  

#exit()

pathwaystr = varlist
newpair=tuple(zip(pathwaystr,MW_g))

annual_dict = dict(zip(varlist, newpair))

for key, val in annual_dict.items():
    print('key:',key)
    print('val=',val)
    if 'TBr2P' in key or 'Prod_Br2' in key:
       annual_dict[key]=(val[0],159.80)

    if 'TCl2P' in key or 'Prod_Cl2' in key:
       annual_dict[key]=(val[0],70.90)

    if 'TI2P' in key or 'Prod_I2' in key:
       annual_dict[key]=(val[0],253.80)

print('--------')
for key, val in annual_dict.items():
    print('key:',key)
    print('val=',val)

# Avogadro's number
avogd = 6.022e23

# unit conversion: from [molecules cm-3 s-1] to [Tg year-1]
for key, value in annual_dict.items():

#  [Tg year-1]    [molecules cm-3 s-1]  [m3]                    [cm3/m3]   [year to seconds]               [MW_g]  [g to Gg]
    ref_tmp = ref_ds[key]     *   ref_metds['Met_AIRVOL'] *  1.0e6   * 365.0 * 24.0 * 3600.0 / avogd * value[1] * 1.0e-9
    dev_tmp = dev_ds[key]     *   dev_metds['Met_AIRVOL'] *  1.0e6   * 365.0 * 24.0 * 3600.0 / avogd * value[1] * 1.0e-9

    ref_tmpsum = 0.0
    dev_tmpsum = 0.0
    for ilat in range(0,len(ref_tmp.lat)):
     for ilon in range(0,len(ref_tmp.lon)):
      itropps = ref_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
      ref_tmpsum=ref_tmpsum + ref_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps)-1)).sum().values

      itropps = dev_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
      dev_tmpsum=dev_tmpsum + dev_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps)-1)).sum().values

    annual_dict[key]=[value[0],value[1],ref_tmpsum,dev_tmpsum]
  

#---------------------

keys=['SpeciesConc_Br','SpeciesConc_Br2','SpeciesConc_BrO','SpeciesConc_BrCl','SpeciesConc_IBr',
      'SpeciesConc_HBr','SpeciesConc_BrNO3','SpeciesConc_HOBr','SpeciesConc_BrNO2',
      'SpeciesConc_BrSALA','SpeciesConc_BrSALC']
       #,'SpeciesConc_Cl',
#     'SpeciesConc_Cl2','SpeciesConc_ClO','SpeciesConc_I','SpeciesConc_I2','SpeciesConc_IO']
MW_O3 = [79.90, 159.80, 79.90, 79.90, 79.90,79.90,79.90,79.90,79.90,79.90,79.90]
         #,35.45,
#         70.90, 35.45, 126.90, 253.80, 126.90]

annual_dict2 = dict(zip(keys, MW_O3))

MW_air = 28.96

lat=ref_metds['Met_AD'].lat
lon=ref_metds['Met_AD'].lon

for key,val in annual_dict2.items():
 
  #  [Gg]    [mol/mol]                    [kg]           [kg -> g]                      [g to Gg]
  ref_tmp = ref_concds[key].values     *   ref_metds['Met_AD'].values *  1.0e3  * val / MW_air   *  1.0e-9
  dev_tmp = dev_concds[key].values     *   dev_metds['Met_AD'].values *  1.0e3  * val / MW_air   *  1.0e-9

  troplev=ref_metds['Met_TropLev'].values
  troplev_dev=dev_metds['Met_TropLev'].values

  ref_tmpsum = 0.0
  dev_tmpsum = 0.0
  print(ref_tmp)
  for ilat in range(0,len(lat)):
   for ilon in range(0,len(lon)):
    #itropps = ref_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
    itropps = troplev[ilat,ilon]
    #ref_tmpsum=ref_tmpsum + ref_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps))).sum().values
    ref_tmpsum=ref_tmpsum + np.sum(ref_tmp[0:int(itropps),ilat,ilon])  

    #itropps = dev_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
    #dev_tmpsum=dev_tmpsum + dev_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps))).sum().values
  
    itropps = troplev_dev[ilat,ilon]
    dev_tmpsum=dev_tmpsum + np.sum(dev_tmp[0:int(itropps),ilat,ilon])  
  annual_dict2[key]=(ref_tmpsum,dev_tmpsum)
  #print(key,':',ref_tmpsum,':',dev_tmpsum)
  
#print('Global Tropospheric O$_3$ burden & ',"{:.1f}".format(O3_burden[0]),' & ',"{:.1f}".format(O3_burden[1]),"\\\\")
print('Global tropospheric species burden (Gg Br/Cl/I)')
for key, value in annual_dict2.items():
    print(key,' & ',"{:.2f}".format(value[0]),' & ',"{:.2f}".format(value[1]),"\\\\")

#-------

plt.figure(figsize = (12,12))

From = ['BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]),
        'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
        'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]),
        'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]),
        'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
        'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1]),
        'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]),
        'BrNO$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][1]),
	'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
	'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
	'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]),
	'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
	'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
	'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
        'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]),
	'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]),
	'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]),
	'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]),
	'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]),
	'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
	'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1]),
	'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
	'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1]),
	'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]),
	'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1]),
	'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
	'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
	'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]),
	'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1]),
	'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
	'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
	'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
        'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1])]

To = ['Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]),
      'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]),
      'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]),
      'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]),
      'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]),
      'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
      'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
      'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
      'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]),
      'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
      'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
      'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
      'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
      'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]),
      'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]),
      'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1]),
      'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1]),
      'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]),
      'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]),
      'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]),
      'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]),
      'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]),
      'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]),
      'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]),
      'BrNO$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][1]),
	'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1])]



Prod=['('+"{:.1E}".format(annual_dict['Prod_TBrPFromBrNO3'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromBrNO3'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrPFromBrO'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromBrO'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrPFromBr2'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromBr2'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrPFromBrCl'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromBrCl'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrPFromHOBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromHOBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrPFromIBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromIBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrPFromHBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromHBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrPFromBrNO2'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrPFromBrNO2'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBr2PFromBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBr2PFromBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBr2PFromHOBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBr2PFromHOBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBr2PFromBrCl'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBr2PFromBrCl'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBr2PFromBrO'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBr2PFromBrO'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrOPFromBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrOPFromBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrOPFromHOBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrOPFromHOBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrOPFromHBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrOPFromHBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrOPFromBrNO3'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrOPFromBrNO3'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THOBrPFromBrNO3'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THOBrPFromBrNO3'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THOBrPFromBrClSS'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THOBrPFromBrClSS'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THOBrPFromBr2SS'][2]+annual_dict['Prod_THOBrPFromBr2wOH'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THOBrPFromBr2SS'][3]+annual_dict['Prod_THOBrPFromBr2wOH'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THOBrPFromBrO'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THOBrPFromBrO'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THOBrPFromO3'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THOBrPFromO3'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrNO3PFromBrO'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrNO3PFromBrO'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TIBrPFromICl'][2]+annual_dict['Prod_TIBrPFromHOI'][2]+annual_dict['Prod_TIBrPFromIONO'][2]+annual_dict['Prod_TIBrPFromIONO2'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TIBrPFromICl'][3]+annual_dict['Prod_TIBrPFromHOI'][3]+annual_dict['Prod_TIBrPFromIONO'][3]+annual_dict['Prod_TIBrPFromIONO2'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TIBrPFromBr2'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TIBrPFromBr2'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrClPFromClNO3'][2]+annual_dict['Prod_TBrClPFromCl2'][2]+annual_dict['Prod_TBrClPFromHOCl'][2]+annual_dict['Prod_TBrClPFromClNO2'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrClPFromClNO3'][3]+annual_dict['Prod_TBrClPFromCl2'][3]+annual_dict['Prod_TBrClPFromHOCl'][3]+annual_dict['Prod_TBrClPFromClNO2'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrClPFromHOBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrClPFromHOBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrClPFromBrO'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrClPFromBrO'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrClPFromBrNO3'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrClPFromBrNO3'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THBrPFromSSBrm'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THBrPFromSSBrm'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THBrPFromBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THBrPFromBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_THBrPFromHOBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_THBrPFromHOBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_TBrNO2PFromBr'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_TBrNO2PFromBr'][3]),
      '('+"{:.1E}".format(annual_dict['Prod_HBrToSS'][2])+')\n'+"{:.1E}".format(annual_dict['Prod_HBrToSS'][3])]


df = pd.DataFrame({ 'source':From,
                   'target':To,
		   'prod':Prod,})
# Define Node Positions
pos = {'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1]):(0.5,0.5+0.5),
        'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]):(0.5-0.3,2+0.3),
        'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]):(0.5,3.5),
        'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]):(2.5,1),
        'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]):(2,2+0.15),
        'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]):(2,3.5+0.15),
        'BrNO$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][1]):(3.5,0.5),
        'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]):(3.5+0.3,2),
        'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]):(3.5,3.5),
	'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1]):(2,0.2)}


# Define Node Colors

NodeColors = {'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1]):[0,0,1],
        'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]):[0.5,0,1],
        'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]):[0,1,1],
        'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]):[0.5,0.5,1],
        'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]):[0.7,0,0],
        'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]):[1,1,0],
        'BrNO$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][1]):[1,0.5,0],
        'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]):[1,0.5,0.5],
        'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]):[0,1,0.5],
	'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1]):[.8,0.2,0.5]}


Labels = {}
i = 0
for a in From:
    Labels[a]=a
    i +=1
Labels[To[-1]]=To[-1]

# Build your graph. Note that we use the DiGraph function to create the graph! This adds arrows
#G=nx.from_pandas_edgelist(df, 'from', 'to','prod', create_using=nx.DiGraph(), edge_attr=True )
G=nx.from_pandas_edgelist(df, create_using=nx.DiGraph(), edge_attr=True )

# Define the colormap and set nodes to circles, but the last one to a triangle
Circles = []
size_Circles = []
Colors_Circles = []
for n in G.nodes:
        Circles.append(n)
        print(n,NodeColors[n])
        Colors_Circles.append(NodeColors[n])
        tmpstr=n.split('\n')
        size_Circles.append(min(((float(tmpstr[2])-0.1)*0.2/(9.0-0.1)+0.1),0.3)*1.0e4)

# By making a white node that is larger, I can make the arrow "start" beyond the node
nodes = nx.draw_networkx_nodes(G, pos, 
                       nodelist = Circles,
                       node_size=size_Circles,
                       node_shape='o',
                       node_color='white',
                       alpha=1)

nodes = nx.draw_networkx_nodes(G, pos, 
                       nodelist = Circles,
                       node_size=size_Circles,
                       node_shape='o',
                       node_color=Colors_Circles,
                       edgecolors='white',
                       alpha=.75)



nx.draw_networkx_labels(G, pos, Labels, font_size=10)

TBrP = []
TBrPlabels = {}
TBrPwidth = []
TBr2P = []
TBr2Plabels = {}
TBr2Pwidth = []

TBrOP = []
TBrOPlabels = {}
TBrOPwidth = []

THOBrP = []
THOBrPlabels = {}
THOBrPwidth = []

TBrNO3P = []
TBrNO3Plabels = {}
TBrNO3Pwidth = []

TIBrP = []
TIBrPlabels = {}
TIBrPwidth = []

TBrClP = []
TBrClPlabels = {}
TBrClPwidth = []

THBrP = []
THBrPlabels = {}
THBrPwidth = []

TBrNO2P = []
TBrNO2Plabels = {}
TBrNO2Pwidth = []

TSSBrP = []
TSSBrPlabels = {}
TSSBrPwidth = []

arc_weight=nx.get_edge_attributes(G,'prod')

for u,v,d in G.edges(data=True):
    print(v)
    if v == 'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

            TBr2Plabels.update({(u,v):d['prod']})
            TBr2P.append((u,v))
            TBr2Pwidth.append(max(math.log10(float(tmpstr[1])),1))
            print('TBr2P',TBr2Plabels)
    elif v == 'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

            TBrPlabels.update({(u,v):d['prod']})
            TBrP.append((u,v))
            TBrPwidth.append(max(math.log10(float(tmpstr[1])),1))
            print('TBrP',TBrPlabels)

    elif v == 'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

            TBrOPlabels.update({(u,v):d['prod']})
            TBrOP.append((u,v))
            TBrOPwidth.append(max(math.log10(float(tmpstr[1])),1))
            print('TBrOP',TBrOPlabels)

    elif v == 'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

            THOBrPlabels.update({(u,v):d['prod']})
            THOBrP.append((u,v))
            THOBrPwidth.append(max(math.log10(float(tmpstr[1])),1))
            print('THOBrP',THOBrPlabels)

    elif v == 'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

           TIBrPlabels.update({(u,v):d['prod']})
           TIBrP.append((u,v))
           TIBrPwidth.append(max(math.log10(float(tmpstr[1])),1))
           print('TIBrP',TIBrPlabels)

    elif v == 'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

           TBrClPlabels.update({(u,v):d['prod']})
           TBrClP.append((u,v))
           TBrClPwidth.append(max(math.log10(float(tmpstr[1])),1))
           print('TBrClP',TBrClPlabels)

    elif v == 'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

           THBrPlabels.update({(u,v):d['prod']})
           THBrP.append((u,v))
           THBrPwidth.append(max(math.log10(float(tmpstr[1])),1))
           print('THBrP',THBrPlabels)

    elif v == 'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

           TSSBrPlabels.update({(u,v):d['prod']})
           TSSBrP.append((u,v))
           TSSBrPwidth.append(max(math.log10(float(tmpstr[1])),1))
           print('TSSBrP',TSSBrPlabels)

    elif v == 'BrNO$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

           TBrNO2Plabels.update({(u,v):d['prod']})
           TBrNO2P.append((u,v))
           TBrNO2Pwidth.append(max(math.log10(float(tmpstr[1])),1))
           print('TBrNO2P',TBrNO2Plabels)

    elif v == 'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1]):
          tmpstr=d['prod'].split('\n')
          if float(tmpstr[1]) > 100:

           TBrNO3Plabels.update({(u,v):d['prod']})
           TBrNO3P.append((u,v))
           TBrNO3Pwidth.append(max(math.log10(float(tmpstr[1])),1))
           print('TBrNO3P',TBrNO3Plabels)

    else:
       print('error exit')
       exit()
# Again by making the node_size larer, I can have the arrows end before they actually hit the node



n = 'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1])

nx.draw_networkx_edges(G, pos,edgelist=TBrP, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=TBrPwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1])

nx.draw_networkx_edges(G, pos,edgelist=TBrOP, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=TBrOPwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1])

nx.draw_networkx_edges(G, pos,edgelist=THOBrP, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=THOBrPwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1])

nx.draw_networkx_edges(G, pos,edgelist=TIBrP, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=TIBrPwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1])

#nx.draw_networkx_edges(G, pos,edgelist=TBrClP, node_size=0.6e4,edge_color=NodeColors[n],
nx.draw_networkx_edges(G, pos,edgelist=TBrClP, node_size=0.6e4,edge_color=[0,0,0],
                               arrowstyle='->',width=TBrClPwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1])

nx.draw_networkx_edges(G, pos,edgelist=TSSBrP, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=TSSBrPwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1])
print('HBr edge color=',NodeColors[n])
print('HBr edges:',THBrP)
print('HBr edges width:',THBrPwidth)
#nx.draw_networkx_edges(G, pos,edgelist=THBrP, node_size=0.6e4,edge_color=NodeColors[n],
nx.draw_networkx_edges(G, pos,edgelist=THBrP, node_size=0.6e4,edge_color=[0,0,0],
                               arrowstyle='->',width=THBrPwidth,arrowsizes=10)#,connectionstyle='arc3,rad=0.2')

n = 'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1])

nx.draw_networkx_edges(G, pos,edgelist=TBrNO3P, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=TBrNO3Pwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'BrNO$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][1])
print('BrNO2 edge color=',NodeColors[n])
print('BrNO2 edges:',TBrNO2P)
print('BrNO2 edges width:',TBrNO2Pwidth)
nx.draw_networkx_edges(G, pos,edgelist=TBrNO2P, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=TBrNO2Pwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1])
nx.draw_networkx_edges(G, pos,edgelist=TBr2P, node_size=0.6e4,edge_color=NodeColors[n],
                               arrowstyle='->',width=TBr2Pwidth,arrowsizes=10,connectionstyle='arc3,rad=0.2')

n = 'Br\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TBrPlabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'BrO\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrO'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TBrOPlabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'BrNO$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO2'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TBrNO2Plabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'BrNO$_{3}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrNO3'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TBrNO3Plabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'HOBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HOBr'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=THOBrPlabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'IBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_IBr'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TIBrPlabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'BrCl\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrCl'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TBrClPlabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'SSBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][0]+annual_dict2['SpeciesConc_BrSALC'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_BrSALA'][1]+annual_dict2['SpeciesConc_BrSALC'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TSSBrPlabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'HBr\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_HBr'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=THBrPlabels,font_color='k',
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

n = 'Br$_{2}$\n'+'('+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][0])+')\n'+"{:.2f}".format(annual_dict2['SpeciesConc_Br2'][1])
nx.draw_networkx_edge_labels(G, pos,edge_labels=TBr2Plabels,font_color='k', 
			    bbox=dict(facecolor=NodeColors[n],edgecolor=[1,1,1], alpha=1), 
			    verticalalignment='center',label_pos=0.70,font_size=6)

plt.xlim(0,4.5)
plt.ylim(0,4)
plt.axis('off')
#plt.show()

plt.savefig('Bry_budget.pdf')





#-----------------------
#keys=['SpeciesConc_O3','SpeciesConc_Br','SpeciesConc_Br2','SpeciesConc_BrO','SpeciesConc_BrCl','SpeciesConc_IBr','SpeciesConc_HBr','SpeciesConc_BrNO3','SpeciesConc_HOBr','SpeciesConc_BrNO2','SpeciesConc_Cl',
#     'SpeciesConc_Cl2','SpeciesConc_ClO','SpeciesConc_I','SpeciesConc_I2','SpeciesConc_IO']
#MW_O3 = [48.0, 79.90, 159.80, 79.90, 79.90, 79.90,79.90,79.90,79.90,79.90,35.45,
#         70.90, 35.45, 126.90, 253.80, 126.90]
#
#annual_dict2 = dict(zip(keys, MW_O3))
#
#MW_air = 28.96
#for key,val in annual_dict2.items():
#  #  [Tg]    [mol/mol]                    [kg]           [kg -> g]                      [g to Tg]
#  #ref_tmp = ref_concds[key]     *   ref_metds['Met_AD'] *  1.0e3  * val / MW_air   *  1.0e-12
#  #dev_tmp = dev_concds[key]     *   dev_metds['Met_AD'] *  1.0e3  * val / MW_air   *  1.0e-12
# 
#  #  [Gg]    [mol/mol]                    [kg]           [kg -> g]                      [g to Gg]
#  ref_tmp = ref_concds[key]     *   ref_metds['Met_AD'] *  1.0e3  * val / MW_air   *  1.0e-9
#  dev_tmp = dev_concds[key]     *   dev_metds['Met_AD'] *  1.0e3  * val / MW_air   *  1.0e-9
# 
#  ref_tmpsum = 0.0
#  dev_tmpsum = 0.0
#  for ilat in range(0,len(ref_tmp.lat)):
#   for ilon in range(0,len(ref_tmp.lon)):
#    itropps = ref_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
#    ref_tmpsum=ref_tmpsum + ref_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps))).sum().values
#  
#    itropps = dev_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
#    dev_tmpsum=dev_tmpsum + dev_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps))).sum().values
#  
#  annual_dict2[key]=(ref_tmpsum,dev_tmpsum)
#  print(key,':',ref_tmpsum,':',dev_tmpsum)
#  
##print('Global Tropospheric O$_3$ burden & ',"{:.1f}".format(O3_burden[0]),' & ',"{:.1f}".format(O3_burden[1]),"\\\\")
#print('Global tropospheric species burden (Gg Br/Cl/I)')
#for key, value in annual_dict2.items():
#    print(key,' & ',"{:.2f}".format(value[0]),' & ',"{:.2f}".format(value[1]),"\\\\")
#
#  
##print('Global Tropospheric O$_3$ burden & ',"{:.1f}".format(O3_burden[0]),' & ',"{:.1f}".format(O3_burden[1]),"\\\\")
#print('#########')
#for key, value in annual_dict.items():
#    print(value[0],' & ',"{:.1f}".format(value[2]),' & ',"{:.1f}".format(value[3]),"\\\\")
#    if key in ['POx_LOx','Prod_Ox','Loss_Ox','PHOx_LHOx','Prod_HOx','Loss_HOx','Prod_TOxPOTHPTW','Prod_TOxLOTHPTW','Prod_THOxPOTHPTW']:
#       print("\\hline")
#
#exit()
#gcplot.compare_single_level(ref_ds, '12.9.3 (Jul 2019)', dev_ds, 'New (Jul 2019)', ilev=ilev, varlist=varlist, pdfname='ProdLoss_12.9.3_vs_new_201907_sfc.pdf')
#gcplot.compare_zonal_mean(ref_ds, '12.9.3 (Jul 2019)', dev_ds, 'New (Jul 2019)', pres_range=[0, 1300], varlist=varlist, pdfname='ProdLoss_12.9.3_vs_new_201907_zonal_mean.pdf')





