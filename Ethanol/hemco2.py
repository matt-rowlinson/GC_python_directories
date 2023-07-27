
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
from datetime import datetime
from netCDF4 import Dataset
import sys
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
plt.style.use("seaborn-darkgrid")

def emission_to_tg(ems, area, unit="kg/m2/s"):
    hold=[] 
    for em in ems:
        x = em * area           ## kg/m2/s to kg/s per gridbox
        x = x * 86400 * 30      ## kg/s/gridbox to kg/gridbox/month
        em =  x * 1e-9          ## kg to Tg
        hold.append( em )
    em = np.sum( np.sum( np.array(hold) , 1 ) , 1 )
    return em

def annual_means( em , times=False):
    try:
        em = np.sum(em.reshape(-1, 12), axis=1)
    except:
        for i in range(1,12):
            try:
                em = np.sum(em[:-i].reshape(-1, 12), axis=1)
                break
            except:
                continue
    if type(times) != bool:
        times = pd.to_datetime(times)
        times = times.year.drop_duplicates().tolist()
        return em, times
    return em


if __name__ == "__main__":
    sectors = ['Total','Anthro', 'Biogenic','Ship','PlantDecay','BioBurn']
    rundir='ceds_only'
    version='13.1.2'

    for sector in sectors:
        NO, times, lat, lon, lev, area, var_list, unit = GC.HEMCO_Diag_read(rundir, version, variable=f'EmisEOH_{sector}', standalone=False)
        if len( NO.shape) == 4:      # If vertical levels then sum over the vertical
            NO = np.nansum( NO, 1 ) 
        NO = emission_to_tg( NO, area )
        NO, time = annual_means( NO, times=times )

        N = NO 
        N_moles = N  / 46.07

        ## plotting
        f, ax = plt.subplots( figsize=(12,4) )
        ax.plot( time[:len(NO)], N_moles, label=f"{sector} Tg N (NO+NH3)" )
        

        plt.ylabel("Tg yr-1")
        plt.legend()
        plt.savefig(f'plots/_{sector}_emissions.png')
        plt.close()
