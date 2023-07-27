#!/usr/bin/env python3
#SBATCH --job-name=inventory_comp_plot
#SBATCH --ntasks=1
#SBATCH --mem=500mb
#SBATCH --partition=interactive
#SBATCH --time=00:07:00
#SBATCH --output=Logs/inventory_comp_plot_%A.py.log
import sys
sys.path.append('/users/mjr583/python_lib/')
import read as R
import RowPy as rp
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

CEDS, Clon, Clat = R.read_em('CEDS', variable='C2H6',  year='2010')
Tz, Tlon, Tlat = R.read_em('tzompa', variable='C2H6')
Ed, Elon, Elat = R.read_em('EDGAR', variable='C2H6')
fltr = Elon>=180.
Elon = np.concatenate((Elon[fltr] - 360.,Elon[~fltr] ))
Ed = np.concatenate((Ed[:,fltr],Ed[:,~fltr]), axis=1)

CEDS_totals = [
        R.regional_totals(CEDS,Clat, Clon,  region='Europe'),
        R.regional_totals(CEDS,Clat, Clon,  region='Asia'),
        R.regional_totals(CEDS,Clat, Clon,  region='North America'),
        R.regional_totals(CEDS,Clat, Clon,  region='South America'),
        R.regional_totals(CEDS,Clat, Clon,  region='Africa'),
        R.regional_totals(CEDS,Clat, Clon,  region='Oceania')
        ]
CEDS_sum = np.sum(CEDS_totals)
CEDS_totals.insert(0, np.nansum(CEDS))
CEDS_totals.insert(-1, np.nansum(CEDS) - CEDS_sum)
#CEDS_totals.insert(1, R.regional_totals(CEDS,Clat,Clon, region='UK') )

Tz_totals = [
        R.regional_totals(Tz,Tlat, Tlon,  region='Europe'),
        R.regional_totals(Tz,Tlat, Tlon,  region='Asia'),
        R.regional_totals(Tz,Tlat, Tlon,  region='North America'),
        R.regional_totals(Tz,Tlat, Tlon,  region='South America'),
        R.regional_totals(Tz,Tlat, Tlon,  region='Africa'),
        R.regional_totals(Tz,Tlat, Tlon,  region='Oceania')
        ]
Tzsum = np.sum(Tz_totals)
Tz_totals.insert(0, np.nansum(Tz))
Tz_totals.insert(-1, np.nansum(Tz) - Tzsum)
#Tz_totals.insert(1, R.regional_totals(Tz,Tlat, Tlon, region='UK') )

Ed_totals = [
        R.regional_totals(Ed,Elat, Elon,  region='Europe'),
        R.regional_totals(Ed,Elat, Elon,  region='Asia'),
        R.regional_totals(Ed,Elat, Elon,  region='North America'),
        R.regional_totals(Ed,Elat, Elon,  region='South America'),
        R.regional_totals(Ed,Elat, Elon,  region='Africa'),
        R.regional_totals(Ed,Elat, Elon,  region='Oceania')
        ]
Edsum = np.sum(Ed_totals)
Ed_totals.insert(0, np.nansum(Ed))
Ed_totals.insert(-1, np.nansum(Ed) - Edsum)
#Ed_totals.insert(1, R.regional_totals(Ed,Elat, Elon, region='UK') )

print(CEDS_totals)
print(Tz_totals)
print(Ed_totals)
NEI = [np.nan, np.nan, np.nan, 0.756*1.25, np.nan ,np.nan ,np.nan ,np.nan]
CEDS_NEI = [CEDS_totals[0]+0.756*1.25, np.nan, np.nan, np.nan, np.nan ,np.nan ,np.nan ,np.nan]

EMEP = [np.nan, 0.0968, np.nan, np.nan, np.nan ,np.nan ,np.nan ,np.nan]


f, ax = plt.subplots(figsize=(10,6))
xlabels = ['Global', 'Europe','Asia','North\nAmerica','South\nAmerica','Africa','Oceania', 'Rest of\nworld']
ax.scatter( xlabels, CEDS_totals, s=100, marker='x', color='#7fc97f',label='CEDS')
ax.scatter( xlabels, Tz_totals, s=100, marker='D', color='#beaed4', label='Tzompa')
ax.scatter( xlabels, Ed_totals, s=100, marker='*', color='#386cb0', label='EDGAR')

ax.scatter( xlabels,NEI, s=100, marker='_', color='grey', label='NEI')
ax.scatter( xlabels,CEDS_NEI, s=100, marker='3', color='grey', label='CEDS + NEI')

#ax.scatter( xlabels,EMEP, s=100, marker='_', color='grey', label='_nolegend_')

ax.set_ylim(bottom=0.)
ax.set_ylabel('Anthropogenic $C_2H_6$ emission (Tg yr$^{-1}$)')
for i in range(len(xlabels)):
    if (i % 2) != 0:
        plt.axvspan(i-.5, i+.5, facecolor='grey', alpha=0.2, zorder=0)
plt.legend(bbox_to_anchor=(.75,.9), fontsize=12, frameon=False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('plots/emissions.png')
plt.close()


f, ax = plt.subplots(figsize=(9,5))
dalsoren = [ CEDS_totals[0]+3.+3.+1.+0.3 ,CEDS_totals[0] , 3.0 ,3.0, 1.0, 0.3 ]
dalsoren_opt = [ 22.3 ,15. , 3.0 ,3.0, 1.0, 0.3 ]

default = [Tz_totals[0], 0., 3.2, 0., 0.]
default.insert(0, np.nansum(default))

dev = [  CEDS_totals[0] , 3., 3.2, 0., 0.]
dev.insert(0, np.nansum(dev))

dev2 = [  CEDS_NEI[0] , 3., 3.2, 0., 0.]
dev2.insert(0, np.nansum(dev2))

optimised=False

xlabels = ['Total','Anthropogenic','Geological','Biomass burning', 'Oceanic','Vegetation']
ax.scatter( xlabels, dalsoren, s=85, marker='v', color='#fdb462',label='Dalsoren et al. (2018)')
if optimised:
    ax.scatter( xlabels, dalsoren_opt, s=85, marker='.', color='#8dd3c7',label='Dalsoren Optimised')
ax.scatter( xlabels, default, s=85, marker='x', color='#bebada', label='Base (Tzompa, no geo)')

ax.scatter( xlabels, dev, s=85, marker='d', color='#fb8072', label='Dev1 (CEDS + Geo)')
ax.scatter( xlabels, dev2, s=85, marker='*', color='#80b1d3', label='Dev2 (CEDS + NEI + Geo)')


ax.set_ylabel('$C_2H_6$ Emission (Tg yr$^{-1}$)')
ax.set_ylim(bottom=0.)
for i in range(len(xlabels)):
    if (i % 2) != 0:
        plt.axvspan(i-.5, i+.5, facecolor='grey', alpha=0.2, zorder=0)
plt.legend(bbox_to_anchor=(0.7,0.87), fontsize=12, frameon=False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
if optimised:
    plt.savefig('plots/sectors_opt.png')
else:
    plt.savefig('plots/sectors.png')

plt.close()




