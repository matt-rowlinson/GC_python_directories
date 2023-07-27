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
import re
from netgraph import Graph
import my_networkx as my_nx
#------------------------------------

refdir='/users/hc2134/scratch/GC/rundirs/merra2_4x5_standard_BryBudgets_Clybudgets_Iybudget_tagged'
devdir='/users/hc2134/scratch/GC/rundirs/merra2_4x5_HBrReversible_Brief_SSHOBrToBr2_SSHOBrToBrCl_SSHOClToCl2_SSHOClToBrCl_BrClToBr2_Cl2ToBrClandICl_Br2ToIBr_HOI_BryBudgets_ClyBudgets_IyBudget_tagged'
#refdir='/users/hc2134/scratch/GC/rundirs/merra2_4x5_HBrReversible_Brief_SSHOBrToBr2_SSHOBrToBrCl_SSHOClToCl2_SSHOClToBrCl_BrClToBr2_Cl2ToBrClandICl_Br2ToIBr_HOI_BryBudgets_ClyBudgets_IyBudget_tagged'

dev_ds = xr.open_dataset(devdir+'/OutputDir/GEOSChem.ProdLoss.20190701_0000z.nc4')
ref_ds = xr.open_dataset(refdir+'/OutputDir/GEOSChem.ProdLoss.20190701_0000z.nc4')


dev_concds = xr.open_dataset(devdir+'/OutputDir_1m/monthly_avg/monthly_avg/GEOSChem.SpeciesConc.20190701_0000z.nc4')
ref_concds = xr.open_dataset(refdir+'/OutputDir_1m/monthly_avg/monthly_avg/GEOSChem.SpeciesConc.20190701_0000z.nc4')


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


ref_ds['Prod_TIClPFromSSCl']=ref_ds['Prod_TIClPFromHOI'] + ref_ds['Prod_TIClPFromIONO'] + ref_ds['Prod_TIClPFromIONO2']
ref_ds.Prod_TIClPFromSSCl.attrs["units"] = "molec cm-3 s-1"

dev_ds['Prod_TIClPFromSSCl']=dev_ds['Prod_TIClPFromHOI'] + dev_ds['Prod_TIClPFromIONO'] + dev_ds['Prod_TIClPFromIONO2']
dev_ds.Prod_TIClPFromSSCl.attrs["units"] = "molec cm-3 s-1"


ref_ds['Prod_TBrClPFromSSCl']=ref_ds['Prod_TBrClPFromHOBr'] 
ref_ds.Prod_TBrClPFromSSCl.attrs["units"] = "molec cm-3 s-1"

dev_ds['Prod_TBrClPFromSSCl']=dev_ds['Prod_TBrClPFromHOBr'] 
dev_ds.Prod_TBrClPFromSSCl.attrs["units"] = "molec cm-3 s-1"


dev_ds['Prod_TIClPFromIOx']=dev_ds['Prod_TIClPFromClOx'] 
dev_ds.Prod_TIClPFromIOx.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TIClPFromIOx']=ref_ds['Prod_TIClPFromClOx'] 
ref_ds.Prod_TIClPFromIOx.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TIClPFromSSI']=ref_ds['Prod_TIClPFromCl2SS'] 
ref_ds.Prod_TIClPFromSSI.attrs["units"] = "molec cm-3 s-1"

dev_ds['Prod_TIClPFromSSI']=dev_ds['Prod_TIClPFromCl2SS'] 
dev_ds.Prod_TIClPFromSSI.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TIBrPFromSSI']=ref_ds['Prod_TIBrPFromBr2'] 
ref_ds.Prod_TIBrPFromSSI.attrs["units"] = "molec cm-3 s-1"


dev_ds['Prod_TIBrPFromSSI']=dev_ds['Prod_TIBrPFromBr2'] 
dev_ds.Prod_TIBrPFromSSI.attrs["units"] = "molec cm-3 s-1"


dev_ds['Prod_TSSIPFromINOx']=dev_ds['Prod_TSSIPFromIONO'] 
dev_ds.Prod_TSSIPFromINOx.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TSSIPFromINOx']=ref_ds['Prod_TSSIPFromIONO'] 
ref_ds.Prod_TSSIPFromINOx.attrs["units"] = "molec cm-3 s-1"

dev_ds['Prod_TSSClPFromBrCl']=dev_ds['Prod_THOBrPFromBrClSS'] 
dev_ds.Prod_TSSClPFromBrCl.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_TSSClPFromBrCl']=ref_ds['Prod_TBrClPFromHOBr'] - ref_ds['Prod_TBrClPFromHOBr']
ref_ds.Prod_TSSClPFromBrCl.attrs["units"] = "molec cm-3 s-1"

ref_ds['Prod_THOClPFromCl2SS']=ref_ds['Prod_THOClPFromCl2SS'] + ref_ds['Prod_THOClPFromCl2wOH']
ref_ds.Prod_THOClPFromCl2SS.attrs["units"] = "molec cm-3 s-1"

dev_ds['Prod_THOClPFromCl2SS']=dev_ds['Prod_THOClPFromCl2SS'] + dev_ds['Prod_THOClPFromCl2wOH']
dev_ds.Prod_THOClPFromCl2SS.attrs["units"] = "molec cm-3 s-1"

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
		     'Prod_TIClPFromSSCl',
		     'Prod_TBrClPFromSSCl',
	   'Prod_TBrClPFromClNO3', 'Prod_TBrClPFromCl2',
           'Prod_TBrClPFromHOCl', 'Prod_TBrClPFromClNO2','Prod_TSSClPFromBrCl',  
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
		     'Prod_TIOxPFromI',
		     'Prod_TIBrPFromICl','Prod_TIBrPFromHOI','Prod_TIBrPFromIONO',
   		     'Prod_TIBrPFromIONO2','Prod_TIBrPFromSSI',
		     'Prod_TIClPFromIONO','Prod_TIClPFromIONO2','Prod_TIClPFromHOI',
		     'Prod_TIClPFromIOx','Prod_TIClPFromSSI','Prod_TSSIPFromHI','Prod_TSSIPFromHOI',
                     'Prod_TSSIPFromIOx','Prod_TSSIPFromINOx','Prod_TSSIPFromIONO2',
     		     'Prod_I','Prod_I2','Prod_IONO2','Prod_INOx',
		     'Prod_HI']


#


#



MW_g_Br =  np.zeros(len(varlistBr))
MW_g_Br[:]=79.90

MW_g_Cl =  np.zeros(len(varlistCl))
MW_g_Cl[:]=35.45

MW_g_I =  np.zeros(len(varlistI))
MW_g_I[:]=126.90


#varlist=varlistBr# + varlistCl + varlistI 
varlist=varlistI #+ varlistI 
#MW_g =np.concatenate(( MW_g_Br, MW_g_Cl, MW_g_I))
#MW_g = MW_g_Br
MW_g = MW_g_I

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

# construct edges uisng production tagging names
edges=[]
G = nx.DiGraph()  # or DiGraph, MultiGraph, MultiDiGraph, etc
#G.add_node("SSCl")
G.add_node("SSI")
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



    if 'Prod_T' in key and 'PFrom' in key:
        tmpstr=re.split('Prod_T|PFrom',key)
        print(tmpstr)
        u=tmpstr[2]
        v=tmpstr[1]


        #if (u in ['HOI','IONO','IONO2']):
        #    continue
        
        #if (u in ['OH','N2O5']):
        #    u = 'SSCl'
        #if (u in ['Cl2SS']):
        #    u = 'Cl2'

        if (u in ['IBrSS']):
            u = 'IBr'

        if (u in ['IClSS']):
            u = 'ICl'

        if (u in ['biINOx']):
            u = 'INOx'

        if (u in ['IONO']):
            u = 'INOx'

        G.add_edge(u, v)
        #if (dev_tmpsum > 100.0): 
        if (dev_tmpsum > 0.0 or ref_tmpsum > 0.0): 
           G[u][v].update({'prod_dev': dev_tmpsum})
           G[u][v].update({'prod_ref': ref_tmpsum})
           G[u][v].update({'label':  '('+np.format_float_scientific(ref_tmpsum,exp_digits=1,precision=1)+')\n'+np.format_float_scientific(dev_tmpsum,exp_digits=1,precision=1)})





print(edges)
#G.add_edges_from(edges)  # using a list of edge tuples
print(G.nodes)
print(G.edges(data=True))
#for u,v,d in G.edges(data=True):
#    print(u,v,d)
#exit()

  


#---------------------
ref_concds['SpeciesConc_ClOx']=ref_concds['SpeciesConc_ClO'] + ref_concds['SpeciesConc_OClO'] + ref_concds['SpeciesConc_ClOO'] + 2.0*ref_concds['SpeciesConc_Cl2O2']
dev_concds['SpeciesConc_ClOx']=dev_concds['SpeciesConc_ClO'] + dev_concds['SpeciesConc_OClO'] + dev_concds['SpeciesConc_ClOO'] + 2.0*dev_concds['SpeciesConc_Cl2O2']

ref_concds['SpeciesConc_IOx']=ref_concds['SpeciesConc_IO'] + ref_concds['SpeciesConc_OIO'] + 2.0*ref_concds['SpeciesConc_I2O2'] + 2.0*ref_concds['SpeciesConc_I2O3'] + 2.0*ref_concds['SpeciesConc_I2O4']  
dev_concds['SpeciesConc_IOx']=dev_concds['SpeciesConc_IO'] + dev_concds['SpeciesConc_OIO'] + 2.0*dev_concds['SpeciesConc_I2O2'] + 2.0*dev_concds['SpeciesConc_I2O3'] + 2.0*dev_concds['SpeciesConc_I2O4']  

ref_concds['SpeciesConc_SSCl']=ref_concds['SpeciesConc_SALACL'] + ref_concds['SpeciesConc_SALCCL'] 
dev_concds['SpeciesConc_SSCl']=dev_concds['SpeciesConc_SALACL'] + dev_concds['SpeciesConc_SALCCL'] 

ref_concds['SpeciesConc_INOx']=ref_concds['SpeciesConc_IONO'] + ref_concds['SpeciesConc_INO'] 
dev_concds['SpeciesConc_INOx']=dev_concds['SpeciesConc_IONO'] + dev_concds['SpeciesConc_INO'] 

ref_concds['SpeciesConc_SSI']=ref_concds['SpeciesConc_ISALA'] + ref_concds['SpeciesConc_ISALC'] 
dev_concds['SpeciesConc_SSI']=dev_concds['SpeciesConc_ISALA'] + dev_concds['SpeciesConc_ISALC'] 
#---------------------

#keys=['SpeciesConc_Br','SpeciesConc_Br2','SpeciesConc_BrO','SpeciesConc_BrCl','SpeciesConc_IBr',
#      'SpeciesConc_HBr','SpeciesConc_BrNO3','SpeciesConc_HOBr','SpeciesConc_BrNO2',
#      'SpeciesConc_BrSALA','SpeciesConc_BrSALC']

#keys=['Cl','Cl2','ClOx',
#      'BrCl','ICl',
#      'HCl','ClNO3','HOCl','ClNO2',
#      'SSCl']
keys=['I','I2','IOx',
      'IBr','ICl',
      'HI','IONO2','HOI','INOx',
      'SSI']

       #,'SpeciesConc_Cl',
#     'SpeciesConc_Cl2','SpeciesConc_ClO','SpeciesConc_I','SpeciesConc_I2','SpeciesConc_IO']
#MW_O3 = [79.90, 159.80, 79.90, 79.90, 79.90,79.90,79.90,79.90,79.90,79.90,79.90]
#MW_O3 = [35.45,       70.90,        35.45,  
#         35.45,       35.45,
#         35.45,       35.45,        35.45,       35.45,
#         35.45]
 
MW_O3 = [126.90,      253.80,        126.90,  
         126.90,      126.90,
         126.90,      126.90,        126.90,       126.90,
         126.90]
         #126.90, 253.80, 126.90]

annual_dict2 = dict(zip(keys, MW_O3))

MW_air = 28.96

lat=ref_metds['Met_AD'].lat
lon=ref_metds['Met_AD'].lon

attrs = {}

for key,val in annual_dict2.items():
 
  #  [Gg]    [mol/mol]                    [kg]           [kg -> g]                      [g to Gg]
  ref_tmp = ref_concds['SpeciesConc_'+key].values     *   ref_metds['Met_AD'].values *  1.0e3  * val / MW_air   *  1.0e-9
  dev_tmp = dev_concds['SpeciesConc_'+key].values     *   dev_metds['Met_AD'].values *  1.0e3  * val / MW_air   *  1.0e-9

  troplev=ref_metds['Met_TropLev'].values
  troplev_dev=dev_metds['Met_TropLev'].values

  ref_tmpsum = 0.0
  dev_tmpsum = 0.0
  print(ref_tmp.shape)
  print(troplev.shape)

  for ilat in range(0,len(lat)):
   for ilon in range(0,len(lon)):
    #itropps = ref_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
    itropps = troplev[0,ilat,ilon]
    #ref_tmpsum=ref_tmpsum + ref_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps))).sum().values
    ref_tmpsum=ref_tmpsum + np.sum(ref_tmp[0,0:int(itropps),ilat,ilon])  

    #itropps = dev_metds['Met_TropLev'].isel(lat=ilat,lon=ilon).values[0]
    #dev_tmpsum=dev_tmpsum + dev_tmp.isel(lat=ilat,lon=ilon,lev=range(0,int(itropps))).sum().values
  
    itropps = troplev_dev[0,ilat,ilon]
    dev_tmpsum=dev_tmpsum + np.sum(dev_tmp[0,0:int(itropps),ilat,ilon])  
  annual_dict2[key]=(ref_tmpsum,dev_tmpsum)
  if key in 'Cl':
      attrs.update({key: {"burden_ref": ref_tmpsum,"burden_dev": dev_tmpsum, 
      "label":key+'\n'+'('+np.format_float_scientific(ref_tmpsum,exp_digits=1,precision=1)+')\n'+np.format_float_scientific(dev_tmpsum,exp_digits=1,precision=1)}})
  else:
      attrs.update({key: {"burden_ref": ref_tmpsum,"burden_dev": dev_tmpsum,
      "label":key+'\n'+'('+"{:.2f}".format(ref_tmpsum)+')\n'+"{:.2f}".format(dev_tmpsum)}})

nx.set_node_attributes(G, attrs)

print(list(G.nodes(data=True)))

#exit()
print('Global tropospheric species burden (Gg Br/Cl/I)')
for key, value in annual_dict2.items():
    print(key,' & ',"{:.2f}".format(value[0]),' & ',"{:.2f}".format(value[1]),"\\\\")

#-------

plt.figure(figsize = (12,12))

# Define Node Positions
#pos = {'ICl':(0.5-0.1,1.3),
#       'ClNO3':(0.3,3.5),
#       'HCl':(2.9,1),
#        'Cl':(2,2+0.3),
#        'ClOx':(2,3.5+0.15),
#        'ClNO2':(3.5,0.7+1),
#        'HOCl':(3.5+0.3,2+0.5),
#        'BrCl':(0.5,2+0.5),
#        'Cl2':(3.5,3.5),
#	'SSCl':(1.5,1)}


pos = {'IBr':(0.5-0.1,1.3),
       'IONO2':(0.3,3.5),
       'HI':(2.9,1),
        'I':(2,2+0.3),
        'IOx':(2,3.5+0.15),
        'INOx':(3.5,0.7+1),
        'HOI':(3.5+0.3,2+0.5),
        'ICl':(0.5,2+0.5),
        'I2':(3.5,3.5),
	'SSI':(1.5,1)}
Labels=nx.get_node_attributes(G,'label')

#--------
# Define the colormap and set nodes to circles, but the last one to a triangle
Circles = []
size_Circles = []
Colors_Circles = []
for n in G.nodes():
    if n not in ['Cl2wOH']:
        Circles.append(n)
        print(n)
        #size_Circles.append(min(((n['burden_dev']-0.1)*0.2/(9.0-0.1)+0.1),0.3)*1.0e4)
        size_Circles.append(min(((G.nodes[n]['burden_dev']-0.1)*0.2/(9.0-0.1)+0.1),0.3)*1.0e4)

nodes = nx.draw_networkx_nodes(G, pos, 
                       nodelist = Circles,
                       node_size=size_Circles,
                       node_shape='o',
                       node_color='white',
                       edgecolors='k',
                       alpha=1)

nx.draw_networkx_labels(G, pos, Labels, font_size=7)

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


curved_edges = []#edge for edge in G.edges()]


for u,v,d in G.edges(data=True):
        print(u)
        print(v)
        print(d)
        if u not in ['Cl2wOH']:

          tmpstr=d['prod_dev']
          if tmpstr > 100.0 :

           curved_edges.append((u,v))
           TIBrPwidth.append(max(math.log10(float(tmpstr)),1))
           print('TIBrP',TIBrPlabels)

    
# Again by making the node_size larer, I can have the arrows end before they actually hit the node

arc_rad=0.25


nx.draw_networkx_edges(G, pos,edgelist=curved_edges, node_size=0.6e4,edge_color='k',
                               arrowstyle='->',width=TIBrPwidth,arrowsizes=10,connectionstyle=f'arc3,rad={arc_rad}')

#-------------------------
edge_weights = nx.get_edge_attributes(G,'label')
curved_edge_labels = {edge: edge_weights[edge] for edge in curved_edges}
#--------------------------
#-------

my_nx.my_draw_networkx_edge_labels(G, pos,font_size=6, edge_labels=curved_edge_labels,rotate=True,rad = arc_rad,label_pos=0.3)
plt.xlim(0,4.5)
plt.ylim(0,4)
plt.axis('off')
#plt.show()

plt.savefig('Iy_budget.pdf')



