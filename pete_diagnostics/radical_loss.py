import matplotlib
matplotlib.use('Agg')
import sys
import copy
import numpy as np
import pandas as pd
import xarray as xr
import ternary as tr
import os
import matplotlib.pyplot as plt
import datetime as datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Patch

#from globals import *

# ~~~~~~~~~~ GC rate calculations
def GCARR(A, B, C, T):
    GCARR =  A * np.exp(C/T) * ((300./T)**B)
    return GCARR

def GCJPLPR(A0, B0, C0, A1, B1, C1, FV, T, NUM_DEN):
    RLOW  = GCARR(A0,B0,C0,T)*NUM_DEN
    RHIGH = GCARR(A1,B1,C1,T)
    
    XYRAT = RLOW/RHIGH
    BLOG  = np.log10(XYRAT)
    FEXP  = 1. / (1. + BLOG * BLOG)
    GCJPLPR = RLOW*FV**FEXP/(1.+XYRAT)
    return GCJPLPR

def GC_HO2HO2(A0, B0, C0, A1, B1, C1, T, H2O, NUM_DEN):
    R0 =  A0 * np.exp(C0/T) * (300./T)**B0
    R1 =  A1 * np.exp(C1/T) * (300./T)**B1

    GC_HO2HO2 = (R0+R1*NUM_DEN)*(1.+1.4E-21 * H2O * np.exp(2200./T))
    return GC_HO2HO2

def GC_TBRANCH(A0, B0, C0, A1, B1, C1, T):
    R0 = A0 * np.exp(C0/T) * ((300./T)**B0)
    R1 = A1 * np.exp(C1/T) * ((300./T)**B1)

    GC_TBRANCH = R0/(1.+R1)
    return GC_TBRANCH
    
def GC_RO2HO2(A0, B0, C0, A1, B1, C1, T):
    R0 =  A0 * np.exp(C0/T) * ((300./T)**B0)
    R1 =  A1 * np.exp(C1/T) * ((300./T)**B1)

    GC_RO2HO2 = R0*(1-np.exp(-0.245E0*R1))
    return GC_RO2HO2
    
def ARRPLUS(A0, B0, C0, D0, E0, T):
    K0 = A0 * (D0+(T*E0))
    ARRPLUS =  K0 * np.exp(-B0/T) * (T/300.)**C0
     
    ARRPLUS[ARRPLUS< 0.0] = 0.0
    return ARRPLUS

def ARSL1K(AREA, RADIUS, GAMMA, T, NUM_DEN):
    MW   = 33
    SRT  = np.sqrt(T)
    SRMW = np.sqrt(MW)
    
    DFKG  = 9.45E+17/NUM_DEN * SRT * np.sqrt(3.472E-2 + 1./ MW)
    ARSL1K = AREA / ( RADIUS/DFKG + 2.749064E-4*SRMW/(GAMMA*SRT))
    return ARSL1K
    
def HETHO2(GAMMA, T, NUM_DEN, AERO, aero_list):
    # Calculate rate for each aerosol species and combine
    for i in range(len(aero_list)):
        if i == 0:
            HETHO2 = ARSL1K(AERO[:,:,:,0,i], AERO[:,:,:,1,i], GAMMA , T, NUM_DEN)
        else:
            HETHO2 += ARSL1K(AERO[:,:,:,0,i], AERO[:,:,:,1,i], GAMMA , T, NUM_DEN)
    return HETHO2

def VPRESH2O(T):
    BOLTZ    = 1.38064852E-23
    CONSVAP  = 6.1078E+3 / ( BOLTZ * 1E+7 )
    CONSEXP  = 17.2693882 * (T - 273.16) / (T - 35.86)
    VPRESH2O = CONSVAP * np.exp(CONSEXP) / T
    return VPRESH2O

def get_aero(path, day, T, NUM_DEN, gamma, latlon=None):
    #Return uptake rate
    
    # ~~~~~~~~~~ Get Aerosol properties
    with xr.open_dataset(f"{path}/OutputDir/GEOSChem.Aerosols.{day}.nc4") as data: 

        data = data.isel({"lev":0,
                  "lat":slice(3,-3),
                  "lon":slice(3,-3)})
        
        if latlon != None:
            data = data.sel({"lat":slice(latlon[0], latlon[2]),
                             "lon":slice(latlon[1], latlon[3])})  
        
        aero_list = ['SSC',
                     'SSA',
                     'OC',
                     'BC',
                     'MDUST7',
                     'MDUST6',
                     'MDUST5',
                     'MDUST4',
                     'MDUST3',
                     'MDUST2',
                     'MDUST1',
                     'SULF']
        
        AERO = []
        for i,x in enumerate(aero_list):
            AERO.append(np.stack([data[f"Chem_AeroArea{x}"].values[:,:,:],   #cm2/cm-3
                                  data[f"Chem_AeroRadi{x}"].values[:,:,:]],  #cm
                                  axis=3))
        AERO = np.stack(AERO, axis=-1)
        HET = HETHO2(float(gamma), T, NUM_DEN, AERO, aero_list)
        return HET
 
# ~~~~~~~~~~ Return model rates
def rate_ratios(path, day, T, NUM_DEN, gamma, latlon=None):
    with xr.open_dataset(f"{path}/OutputDir/GEOSChem.ConcAfterChem.{day}.nc4") as data: 
        
        data = data.isel({"lev":0,
                          "lat":slice(3,-3),
                          "lon":slice(3,-3)})
                          
        if latlon != None:
            data = data.sel({"lat":slice(latlon[0], latlon[2]),
                             "lon":slice(latlon[1], latlon[3])})                
                      
        HO2_CONC = data["concAfterChem_HO2"].values[:,:,:]
        OH_CONC  = data["concAfterChem_OH"].values[:,:,:]
        NO2_CONC = data["concAfterChem_NO2"].values[:,:,:]
        H2O = data["concAfterChem_H2O"].values[:,:,:]
        
        #NO2 + OH rate
        NO2 = GCJPLPR(1.80E-30, 3.0E+00, 0.0, 2.8E-11, 0.0, 0.0, 0.6, T, NUM_DEN) * NO2_CONC * OH_CONC

        ROX = 0
        #ROx + ROx
        with open("/users/mjr583/GC/pete_diagnostics/All_reactions_radical_loss.txt") as file:
            LINES = file.readlines()
            for i,L in enumerate(LINES):
                REACTS = L[:L.find("=")].split()
                A = data[f"concAfterChem_{REACTS[0]}"].values[:,:,:]
                try:
                    B = data[f"concAfterChem_{REACTS[2]}"].values[:,:,:]
                except:
                    continue
                #if np.sum(A) == 0.:
                #    print('A', REACTS)
                #elif np.sum(B) == 0.:
                #    print('B', REACTS)
                if "NITs" in REACTS:
                    print( REACTS, A.sum(), B.sum() )
                    sys.exit()

                REACTION = L[L.find(":")+1: L.find("!")].strip()
                #Replace the brackets with comas for delimiting
                REACTION = REACTION.replace(")", ",").replace("(", ",").split(",")
                #REACTION = list(filter(None, REACTION)) # Drop empty list elements 
                REACTION = [ s.replace("d","E") for s in REACTION ] # Replace "d" double precision with python readable "E"

                #print(REACTION)
                if "GCARR" in REACTION[0]:
                    if "_abc" in REACTION[0]:
                        RATE = GCARR(float(REACTION[1]), float(REACTION[2]), float(REACTION[3]), T) * A * B
                    elif "_ac" in REACTION[0]:
                        RATE = GCARR(float(REACTION[1]), 0., float(REACTION[2]), T) * A * B
                elif "GCJPLPR" in REACTION[0]:
                    if "_aba" in REACTION[0]:
                        RATE = GCJPLPR(float(REACTION[1]), float(REACTION[2]), 0.0, float(REACTION[3]), 0.0, 0.0, float(REACTION[4]), T, NUM_DEN) * A * B
                elif "GC_HO2HO2" in REACTION[0]:
                    RATE = GC_HO2HO2(float(REACTION[1]), 0., float(REACTION[2]), float(REACTION[3]), 0., float(REACTION[4]) , T, H2O, NUM_DEN) * A * B
                elif REACTION[0] == "GC_TBRANCH_1_acac":
                    RATE = GC_TBRANCH(float(REACTION[1]), 0., float(REACTION[2]), float(REACTION[3]), 0., float(REACTION[4]), T) * A * B
                elif REACTION[0] == "GC_RO2HO2_aca":
                    RATE = GC_RO2HO2(float(REACTION[1]), 0., float(REACTION[2]), float(REACTION[3]), 0., 0., T) * A * B

                elif REACTION[0] == "ARRPLUS_abEe":
                    RATE = ARRPLUS(float(REACTION[1]), float(REACTION[2]), 0., float(REACTION[3]), float(REACTION[4]), T) * A * B

                elif isfloat(REACTION[0]) == True:
                    RATE = float(REACTION[0])
                else:
                    print(REACTION[0], "MISSING")
                RADICALS_LOST = L[L.find("!")+17:].strip().rstrip(")") 
                
                RATE = RATE * -float(RADICALS_LOST)
                if i == 0:
                    ROX = copy.copy(RATE)
                else:
                    ROX += RATE

        #HO2 uptake aerosol
        HET = get_aero(path, day, T, NUM_DEN, gamma, latlon=latlon) * HO2_CONC 
           
        #HET[:,:,:] = 0. #!!!!!!!!!!!! DELETE THIS WHEN NOT USING FOR PLOT FOR MAT
        #Total Radical loss
        TOTAL = NO2+ROX+HET
    
    return NO2, ROX, HET, TOTAL
    
def LON_NAX_ZENITH(data, SUN):
    NEW_DATA = []
    LON_SUN_TIME = np.argmax(SUN[:,0,:], axis=0)
    for i in range(data.shape[2]):
        NEW_DATA.append(data[LON_SUN_TIME[i],:,i])
    NEW_DATA = np.array(NEW_DATA).T
    return NEW_DATA
    
def add_box(ax, bounds):
    ax.plot([bounds[1],bounds[3]], [bounds[0],bounds[0]], color="k", linewidth=0.75)
    ax.plot([bounds[1],bounds[3]], [bounds[2],bounds[2]], color="k", linewidth=0.75)
    ax.plot([bounds[1],bounds[1]], [bounds[0],bounds[2]], color="k", linewidth=0.75)
    ax.plot([bounds[3],bounds[3]], [bounds[0],bounds[2]], color="k", linewidth=0.75)

def calculate_day(path, D, gamma, latlon):
    #Day
    day = f"201601{D:02}_0000z"
    print(day)
    # ~~~~~~~~~~ Get Met variables and dimension indexes
    with xr.open_dataset(f"{path}/OutputDir/GEOSChem.StateMet.{day}.nc4") as data: 

        data = data.isel({"lev":0,
                          "lat":slice(3, -3),
                          "lon":slice(3, -3)})
        
        if latlon != None:
            data = data.sel({"lat":slice(latlon[0], latlon[2]),
                             "lon":slice(latlon[1], latlon[3])})
        
        LAT = data.lat.values
        LON = data.lon.values
        
        LON = np.append(LON-((LON[1]-LON[0])/2), LON[-1]+((LON[1]-LON[0])/2))
        LAT = np.append(LAT-((LAT[1]-LAT[0])/2), LAT[-1]+((LAT[1]-LAT[0])/2))
          
        SUN = data["Met_SUNCOS"].values[:,:,:]
        T = data["Met_T"].values[:,:,:]

        NUM_DEN = (data["Met_AIRDEN"].values[:,:,:]/1E6)*1E3 #g/cm-3
        NUM_DEN = (NUM_DEN / 28.97) * 6.0221409E+23 #Molc Air/cm-3
     #take mean value for calculation as there is minimal change
     
    NO2, ROX, HET, TOTAL = rate_ratios(path, day, T, NUM_DEN, gamma, latlon=latlon)
        
    TOTAL = LON_NAX_ZENITH(TOTAL, SUN)
    NO2 = LON_NAX_ZENITH(NO2, SUN)
    ROX = LON_NAX_ZENITH(ROX, SUN)
    HET = LON_NAX_ZENITH(HET, SUN)
    data = np.array([NO2, ROX , HET, TOTAL])
    return data, LAT, LON

# EVERYTHING FUNCTION ABOVE THIS POINT DEALS WITH WORKING UP THE RATES. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Plot ~~~~~~~~~~~~~~~~~~~~
def RGB_plot_year(ax, path, year, n_bins, name="", gamma=0.2, latlon=None):
    
    #Change the year label to match the emission year rather than the met year if year is 2017.
    if year == "2017":
        year = "2014"
        
    data_stack = []
    for D in range(1, 4):
        print(year, D)
        data, LAT, LON = calculate_day(path, D, gamma, latlon = latlon)
        data_stack.append(data)
    data = np.stack(data_stack, axis=-1).mean(axis=-1)
    #Divide by total
    data = data/data[3,:,:]    
    RGB = data[:3,:,:] 
    
    proj = ccrs.PlateCarree()

    imshow_extent = (LON[0], LON[-1], LAT[0], LAT[-1])
    
    if n_bins == 2:
        # 2 bins is a special case. If no species is over 50% then all species get rounded to 0 and thus black.
        # For 2 bins just round the closet corner
        RGB = np.floor(RGB/RGB.max(axis=0))
        
        #Print number of boxes in each regime
        print(f"VOC limited {(np.count_nonzero(RGB[0,:,:])/RGB[0,:,:].size):.3f}")
        print(f"NOx limited {(np.count_nonzero(RGB[1,:,:])/RGB[1,:,:].size):.3f}")
        print(f"Aerosol inhibited {(np.count_nonzero(RGB[2,:,:])/RGB[2,:,:].size):.3f}")
        
    else:
        # Used with crop normalize to move regieme in to the middle of colors  
        RGB = np.floor(RGB*n_bins)
        RGB[RGB==n_bins] = n_bins-1
        RGB = RGB*(1/(n_bins-1))
        
        # Crop color normalize
        # https://stackoverflow.com/questions/54056392/maxwell-color-triangle-with-matplotlib
        RGB = RGB*3.
        RGB[RGB>1.] = 1.
        
    # Reshape for imshow (Lat, LON, RGB)
    RGB = RGB.transpose(1,2,0)

    ax.imshow(RGB,origin="lower", extent=imshow_extent, transform=proj)
    
    ax.coastlines(resolution="10m", linewidth=0.6)
    
    ax.add_feature(cfeature.NaturalEarthFeature("cultural", "admin_0_countries", "10m"),
                                                facecolor="none",
                                                edgecolor='black',
                                                linewidth=0.6)
                                                
    ax.add_feature(cfeature.NaturalEarthFeature("physical", "lakes", "10m"),
                                                facecolor="none",
                                                edgecolor='black',
                                                linewidth=0.6)                                                   
           
     
    #ax.set_xlim(left=LON[0], right=LON[-1])
    #ax.set_ylim(bottom=LAT[0], top=LAT[-1])
    label_loc_lon = LON[0]
    #12% from the top
    label_loc_lat = LAT[-1]-((LAT[-1]-LAT[0])*0.12)
    if name == "":
        ax.text(label_loc_lon, label_loc_lat, f"{year}", backgroundcolor="white", fontsize=22, bbox=dict(boxstyle='square,pad=0', fc='white', ec='none', alpha=0.8))
    else:
        ax.text(label_loc_lon, label_loc_lat, f"{year} ${name}$", backgroundcolor="white", fontsize=22, bbox=dict(boxstyle='square,pad=0', fc='white', ec='none', alpha=0.8))

def hex_plot_CB(ax, n_bins):
    #n_bins must be an odd number 3 or greater
    
    scale = n_bins-1
    tax = tr.TernaryAxesSubplot(ax=ax,scale=scale)

    tax.boundary(linewidth=0.75)
#    tax.gridlines(multiple=1, color="blue")

    bins = np.array(range(n_bins))
    
    data = {}
    for i in bins:
        I = i*(1/(n_bins-1))
        I = I*3
        if I > 1.:
            I = 1.
        for j in bins:
            J =  j*(1/(n_bins-1))
            J = J*3
            if J > 1.:
                J = 1.
            for k in bins:
                K = k*(1/(n_bins-1))
                K = K*3
                if K > 1.:
                    K = 1.
                if (i + j + k)!=scale:
                    continue
                data[(i,j,k)] = (I,K,J,1) # normalise
                
        

    tax.heatmap(data, scale=scale, style="hexagonal", use_rgba=True, colorbar=False) 


    tax.bottom_axis_label(r"$\mathsf{NO_2\ +\ OH\ \rightarrow}$", fontsize=12,color="k", offset=0.03)
    tax.right_axis_label(r"$\mathsf{\leftarrow HO_2\ uptake}$", fontsize=12,  color="k", offset=0.075)
    tax.left_axis_label(r"$\mathsf{\leftarrow RO_x\ +\ RO_x}$", fontsize=12,  color="k", offset=0.075)
     
    #tax.ticks(axis='lbr', linewidth=1, multiple=0.25, tick_formats="%.1f", offset=0.02)
    tax.clear_matplotlib_ticks()

    #ax.set_title(name, fontsize="xx-large")

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)   

def RGB_plot(latlon=None, name="", n_bins=7):
    proj = ccrs.PlateCarree()
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3,ncols=1,figsize=(18, 8), subplot_kw=dict(projection=proj))

    RGB_plot_year(ax1, f"/mnt/lustre/users/mjr583/GC/13.1.2_diagnostics/rundirs/TEST_100-scale/", "Jan", n_bins, latlon=latlon)
    #RGB_plot_year(ax2, f"/mnt/lustre/users/mjr583/GC/13.1.2_diagnostics/rundirs/test/", "2019b", n_bins, latlon=latlon)
    #RGB_plot_year(ax3, f"/mnt/lustre/users/mjr583/GC/13.1.2_diagnostics/rundirs/test/", "2019c", n_bins, latlon=latlon)
    
    fig.subplots_adjust(wspace=0.0, hspace=0.04)
    """
    if latlon==None:
        add_box(ax1, (30, 129, 46, 147))
        add_box(ax1, (11.5, 68, 32, 94))
        add_box(ax1, (30, -125, 50, -66))
        add_box(ax1, (36.5, -13, 58.5, 28))
        add_box(ax1, (20, 100, 43, 125))
    """
    pos = ax3.get_position().bounds
    cax = fig.add_axes([pos[0], pos[1]+0.01, (0.20*8)/18, 0.20])
    hex_plot_CB(cax, n_bins)
    cax.patch.set_alpha(0.2)
    if name != "":
        name = "_"+name    
    plt.savefig(f"Jscale_plots/RGB_{name}_noon.png".replace(" ","_"), dpi=300, bbox_inches="tight")
    plt.close()

def RGB_4_plot_dif(year, run1, run2, run3, name, gamma=0.2):
    proj = ccrs.PlateCarree()
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4,ncols=1,figsize=(18, 10.66), subplot_kw=dict(projection=proj))
    RGB_plot_year(ax1, f"../{year}/geosfp_0.5x0.625_tropchem",                    year, name="\mathsf{base}",  n_bins=7)
    RGB_plot_year(ax2, f"../{year}_Experiments/geosfp_0.5x0.625_tropchem_{run1}", year, name=RUN_DICT[run1],   n_bins=7, gamma=gamma)
    RGB_plot_year(ax3, f"../{year}_Experiments/geosfp_0.5x0.625_tropchem_{run2}", year, name=RUN_DICT[run2],   n_bins=7, gamma=gamma)
    RGB_plot_year(ax4, f"../{year}_Experiments/geosfp_0.5x0.625_tropchem_{run3}", year, name=RUN_DICT[run3],   n_bins=7, gamma=gamma)
    
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    pos = ax4.get_position().bounds
    cax = fig.add_axes([pos[0], pos[1]+0.01, ((0.20*10.66)/18)*0.75, 0.20*0.75])
    hex_plot_CB(cax, n_bins=7)
    cax.patch.set_alpha(0.2)
    
    plt.savefig(f"Jscale_plots/RGB_{year}_{name}.png".replace(" ","_"), dpi=300, bbox_inches="tight")
    plt.close()

def regime_plot(latlon=None, name=""):
    proj = ccrs.PlateCarree()
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3,ncols=1,figsize=(18, 8), subplot_kw=dict(projection=proj))

    RGB_plot_year(ax1, f"/mnt/lustre/users/mjr583/GC/13.1.2_diagnostics/rundirs/TEST_1-scale", "2016", n_bins=2, latlon=latlon)
    #RGB_plot_year(ax2, f"../1970/geosfp_0.5x0.625_tropchem", "1970", n_bins=2, latlon=latlon)
    #RGB_plot_year(ax3, f"../2017/geosfp_0.5x0.625_tropchem", "2014", n_bins=2, latlon=latlon)
    
    fig.subplots_adjust(wspace=0.0, hspace=0.04)
    
    legend_elements = [Patch(facecolor=(1, 0, 0), label='$\mathsf{VOC\ limited}$'),
                       Patch(facecolor=(0, 1, 0), label='$\mathsf{NO_x\ limited}$'),
                       Patch(facecolor=(0, 0, 1), label='$\mathsf{Aerosol\ inhibited}$'),
                       ]
    
    plt.legend(handles=legend_elements, loc="lower left",fontsize="x-large")
    if name != "":
        name = "_"+name
    plt.savefig(f"Jscale_plots/Regime{name}_noon.png".replace(" ","_"), dpi=300, bbox_inches="tight")
    plt.close()
    
def regime_plot_dif(year, runs, name, gamma_list=None):
    proj = ccrs.PlateCarree()
    n_runs = len(runs)+1
    if gamma_list == None:
        gamma_list = [0.2 for r in runs]
        
    fig, axes = plt.subplots(nrows=n_runs, ncols=1, figsize=(18, 2.65*n_runs), squeeze=True, subplot_kw=dict(projection=proj))
    
    RGB_plot_year(axes[0], f"../{year}/geosfp_0.5x0.625_tropchem", year, name="\mathsf{base}", n_bins=2)
    for i, run in enumerate(runs):
        RGB_plot_year(axes[i+1], f"../{year}_Experiments/geosfp_0.5x0.625_tropchem_{run}", year, name=RUN_DICT[run],   n_bins=2, gamma=gamma_list[i])
    
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    
    legend_elements = [Patch(facecolor=(1, 0, 0), label='$\mathsf{VOC\ limited}$'),
                       Patch(facecolor=(0, 1, 0), label='$\mathsf{NO_x\ limited}$'),
                       Patch(facecolor=(0, 0, 1), label='$\mathsf{Aerosol\ inhibited}$'),]
    
    plt.legend(handles=legend_elements, loc="lower left",fontsize="x-large")
    
    plt.savefig(f"Jscale_plots/Regime_{year}{name}.png".replace(" ","_"), dpi=300, bbox_inches="tight")
    plt.close() 

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

if __name__ == "__main__":
    # RGB_plot(name="no_HET")
    RGB_plot(n_bins=9)
    
    #x = ARRPLUS(3.94E-12, 0., 0., -0.2038, 9.0435E-4)
    #print(x)
    #RGB_plot((20,   20,    50, 150),  "Asia")
    #RGB_plot((30,  -50, 58.75,  60),  "Europe")
    #RGB_plot((25, -175, 58.75, -40),  "America")

    # regime_plot(name="no_HET")
    #regime_plot()

    #regime_plot((20,   20,    50, 150),  "Asia")
    #regime_plot((30,  -50, 58.75,  60),  "Europe")
    #regime_plot((25, -175, 58.75, -40),  "America")

    # regime_plot_dif("1970",  ["no_ships"], "no_ships")
    # RGB_4_plot_dif("1750", "no_BB", "no_DST", "no_SS", "natural")
    # regime_plot_dif("2017",  ["GAMMA_01"], "gamma", gamma_list=[0.1])
    # regime_plot_dif("2017",  ["GAMMA_01", "GAMMA_00"], "gamma", gamma_list=[0.1, 1E-20])


