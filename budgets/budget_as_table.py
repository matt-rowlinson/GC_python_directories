import matplotlib.pyplot as plt
import xarray as xr
import os
import sys
import numpy as np
import pandas as pd
#--------------------------------
pos = {
       'pNO3'               : (-7,0),
       'NO'                 : (26,-20),
       'NO2'                : (-8,-31),
       'N2O5'               : (-16,23),
       'HNO3'               : (28,17),
       'PAN'                : (-28, 26),
       'HONO'               : (4,11),
       'Halogen nitrates'   : (16,30),
       'NO3'                : (-28,-25)}

options = {
            "node_size" : 7e2,
            "node_color": "white",
            "edgecolors": "white",
            "alpha"     :0.9,
            "node_shape":'o'}
        
Labels = {'NO'                  :'NO',
          'HONO'                :'HONO',
          'PAN'                 :'PAN',
          'pNO3'                :'pNO$_3^-$',
          'N2O5'                :'N$_2$O$_5$',
          'NO2'                 :'NO$_2$',
          'HNO3'                :'HNO$_3$',
          'Halogen nitrates'    :'(Br,Cl,I)NO$_3$',
          'NO3'                 :'NO$_3$'
          }

e_list = {  ('NO', 'NO2') : {'O3 + NO  --> NO2':    ['NG001'],
                             'HO2 + NO --> NO2':    ['NG003'],
                             'MO2 + NO --> NO2':    ['NG004'],
                             'NO + NO3 --> NO2':    ['NG007'],
                             'Total':               ['NG001','NG003','NG004','NG007']},
            ('NO2', 'NO') : {'NO2 + NO3 --> NO':    ['NG011'],
                             'NO2 + hv  --> NO':    ['OP001'],
                             'Total':               ['NG011','OP001']},
            ('HONO', 'NO2') : { 'HONO + OH --> ':   ['NG005']},
            ('NO3', 'NO2') : {'HO2 + NO3   --> NO2':['NG006'],
                              'NO + NO3    --> NO2':['NG007'],
                              'OH + NO3    --> NO2':['NG008'],
                              'NO3 + hv    --> NO2':['OP004'],                
                              'Total' :             ['NG006','NG007','NG008','OP004']},
            ('NO2','HNO3') : {'NO2 + OH + M --> HNO3 + M' : ['NH001']},
            ('NO2', 'N2O5') : {'NO2 + NO3 --> N2O5':['NG009']},
            ('NO3', 'N2O5') : {'NO2 + NO3 --> N2O5':['NG009']},
            ('N2O5', 'NO2') : {'N2O5 + M  --> NO2 + NO3':['NG010']},
            ('N2O5', 'NO3') : {'N2O5 + hv --> NO2 + NO3':['OP006']},
            ('NO3', 'NO')   : {'NO2 + NO3 --> NO'  :['NG011'],
                               'NO3 + hv  --> NO'  :['OP005'],
                               'Total' : ['NG011','OP005']},
            ('NO2','PAN')   : {'MO2 + NO2 --> PANs':['NG012']},
            ('PAN','NO2')   : {'PANs      --> NO2' :['NG013'],
                               'PAN + hv  --> NO2' :['OP007'],
                               'Total' : ['NG013','OP007']},
            ('NO', 'HONO')  : {'NO + OH --> HONO'            :['NG014']},
            ('HNO3', 'NO3') : {'HNO3 + OH --> H2O + NO3'     :['NG015']},
            ('pNO3','HONO') : {'NIT --> HONO'       :['NP001','NP002','NP003']},
            ('pNO3','NO2')  : {'NIT --> NO2'        :['NP004','NP005','NP006']},
            ('HNO3','NO2')  : {'HNO3 + hv --> OH + NO2'    :['OP002']},
            ('HONO','NO')   : {'HONO + hv --> OH + NO'      :['OP003']},
            ('PAN','NO3')   : {'PAN + hv --> NO3'            :['OP007']},
            ('NO2','HONO')  : {'NO2 --> 0.5 HNO2'            :['HH001']},
            ('NO2','HNO3')  : {'NO2 --> 0.5 HNO3'            :['HH001']},
            ('NO3','HNO3')  : {'NO3 --> HNO3'           :['HH002']},            
            ('N2O5','HNO3') : {'N2O5 + H2O --> 2HNO3'           :['HH003']},
            ('Halogen nitrates','HNO3') : {'Cl/Br/INO3 + H2O/HCl --> HNO3':['HH004','HH005','HH006']},
            ('Halogen nitrates','HONO') : {'Cl/Br/INO2 + H2O/HCl --> HONO':['HN001','HN002']}
            }



#--------------------------------
def find_file_list(path, substrs):
    file_list =[]
    for root, directory, files in os.walk(path):
        for f in files:
            for s in substrs:
                if s in f:
                    file_list.append(os.path.join(root, f))
    file_list.sort()
    return file_list

def get_data_as_xr(rundir, dates='2017', variable=False):
    path=f'/mnt/lustre/users/mjr583/GC/14.1.0/rundirs/{rundir}'
    ds = xr.open_mfdataset( find_file_list(path,dates), combine='by_coords' )
    return ds

def label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                MW, tracers):
    avo=6.022e23
    MW_air = 28.96
    tlist = [ "Prod_" + s for s in tracers]
    
    ds = ref_spc[tlist].to_array().sum("variable")
    ds = ds * ref_met['Met_AIRVOL'] *  1.0e6   * 365.0 * 24.0 * 3600.0 / avo * MW * 1.0e-12
    #ds = ds * ref_met['Met_AD'].values * 1.0e3 / MW_air * MW * 1.0e-9
    Trop_Mask = (ref_met['Met_PMID'] > ref_met['Met_TropP'])
    ds1 = ds.where(Trop_Mask)
    ref_value = ds1.sum().values
    
    ds = dev_spc[tlist].to_array().sum("variable")
    #ds = ds * dev_met['Met_AD'].values * 1.0e3 / MW_air * MW * 1.0e-9
    ds = ds * dev_met['Met_AIRVOL'] *  1.0e6   * 365.0 * 24.0 * 3600.0 / avo * MW * 1.0e-12
    ds1 = ds.where(Trop_Mask)
    dev_value = ds1.sum().values

    return ref_value, dev_value

####------------------------
def main():
     
    ## Read GEOS-Chem tagged run output
    ref_met = get_data_as_xr('new_base_4x5_tracers',dates=["StateMet.2019"]).mean(dim='time')
    ref_spc = get_data_as_xr('new_base_4x5_tracers',dates=["SpeciesConc.2019"]).mean(dim='time')
    ref_pro = get_data_as_xr('new_base_4x5_tracers',dates=["ProdLoss.20190810_0200"]).mean(dim='time')


    dev_met = get_data_as_xr('langmuir_4x5_tracers',dates=["StateMet.2019"]).mean(dim='time')
    dev_spc = get_data_as_xr('langmuir_4x5_tracers',dates=["SpeciesConc.2019"]).mean(dim='time')
    dev_pro = get_data_as_xr('langmuir_4x5_tracers',dates=["ProdLoss.20190810_0200"]).mean(dim='time')

    # Label each pathway one by one using the dictionary of pathways and tracers
    first = True
    for edge, tracers in e_list.items():
        for rxn, tracer in tracers.items():
            if edge[1]=="N2O5":
                MW=14.007*2
            else:
                MW=14.007
            ref_value, dev_value = label_rxns( ref_met, ref_pro,
                                               dev_met, dev_pro,
                                               MW, tracer)
            if first:
                df = pd.DataFrame({"Pathway":f"{edge[0]} --> {edge[1]}",
                                   "Reaction": rxn,
                                   "Base (Tg N)": '%.2f' %ref_value,
                                   "NIThv On (Tg N)": '%.2f' %dev_value}, index=[' '])
                first = False
            else:
                a_df = pd.DataFrame({"Pathway":f"{edge[0]} --> {edge[1]}",
                                     "Reaction":rxn,
                                     "Base (Tg N)": '%.2f' %ref_value, 
                                     "NIThv On (Tg N)": '%.2f' %dev_value},  index=[' '])
                df=pd.concat([df,a_df])
    print( df )
    
    e_df = pd.DataFrame({"Pathway":" --> NOx",
                         "Reaction":'Emission',
                         "Base (Tg N)": '%.2f' %44.0, 
                         "NIThv On (Tg N)": '%.2f' %44},  index=[' '])
    f_df = pd.DataFrame({"Pathway":" --> HNO3",
                         "Reaction":'Emission',
                         "Base (Tg N)": '%.2f' %4.0, 
                         "NIThv On (Tg N)": '%.2f' %4.0},  index=[' '])

    g_df = pd.DataFrame({"Pathway":"HNO3 --> ",
                         "Reaction":'Deposition',
                         "Base (Tg N)": '%.2f' %24.0, 
                         "NIThv On (Tg N)": '%.2f' %24.0},  index=[' '])
    h_df = pd.DataFrame({"Pathway":"PAN --> ",
                         "Reaction":'Deposition',
                         "Base (Tg N)": '%.2f' %2.0, 
                         "NIThv On (Tg N)": '%.2f' %2.0},  index=[' '])
    i_df = pd.DataFrame({"Pathway":"N2O5 --> ",
                         "Reaction":'Deposition',
                         "Base (Tg N)": '%.2f' %0.21, 
                         "NIThv On (Tg N)": '%.2f' %0.25},  index=[' '])
    j_df = pd.DataFrame({"Pathway":"NO2 --> ",
                         "Reaction":'Deposition',
                         "Base (Tg N)": '%.2f' %3.0, 
                         "NIThv On (Tg N)": '%.2f' %3.},  index=[' '])
    k_df = pd.DataFrame({"Pathway":"pNO3 --> ",
                         "Reaction":'Deposition',
                         "Base (Tg N)": '%.2f' %7.0, 
                         "NIThv On (Tg N)": '%.2f' %4.0},  index=[' '])
    l_df = pd.DataFrame({"Pathway":"HONO --> ",
                         "Reaction":'Deposition',
                         "Base (Tg N)": '%.2f' %0.17, 
                         "NIThv On (Tg N)": '%.2f' %0.364},  index=[' '])


    df = pd.concat([ df,e_df,f_df,g_df,h_df,i_df,j_df,k_df,l_df])

    print( df )
    fig, ax =plt.subplots(figsize=(9,5))
    ax.axis('tight')
    ax.axis('off')
    rcolors = plt.cm.BuPu(np.full(len(df.index), 0.1))
    ccolors = plt.cm.BuPu(np.full(len(df.columns), 0.1))

    the_table = ax.table(cellText=df.values,
            rowLabels=df.index,
            rowColours=rcolors,
            rowLoc='right',
            colColours=ccolors,  
            colLabels=df.columns,
            loc='center')
    plt.tight_layout
    plt.savefig('plots/Budget_table_2.png', dpi=300, bbox_inches='tight')
    plt.close()


    sys.exit()


   



    avo=6.022e23
    hemco = get_data_as_xr('langmuir_4x5_tracers',dates=["HEMCO_diag"]).mean(dim='time')
    ref_prodl = get_data_as_xr('langmuir_4x5_tracers',dates=["ProdLoss.2019"]).mean(dim='time')
    ref_drydp = get_data_as_xr('langmuir_4x5_tracers',dates=["DryDep.2019"]).mean(dim='time')

    ## Read emissions totals from HEMCO
    NO_em  = ( hemco['EmisNO_Total'] * hemco['AREA'] * 365. * 24. * 3600 * 1e-9 ).sum().values * (14.007/46.0055)
    NO2_em = ( hemco['EmisNO2_Total'] * hemco['AREA'] * 365. * 24. * 3600 * 1e-9 ).sum().values * (14.007/30.01)
    NOx_em = np.round(NO_em + NO2_em)
    HNO3_em = ( hemco['EmisHNO3_Ship'] * hemco['AREA'] * 365. * 24. * 3600 * 1e-9 ).sum().values.round() * (14.007/63.01)
    
    ## Read dry dep total fro DryDep
    NO2_dd  = ( ref_drydp['DryDep_NO2']  * hemco['AREA'] * 1.0e4 * 365. * 24. * 3600 / avo * 14.007 * 1e-12).sum().values 
    HNO3_dd = ( ref_drydp['DryDep_HNO3'] * hemco['AREA'] * 1.0e4 * 365. * 24. * 3600 / avo * 14.007 * 1e-12).sum().values
    PAN_dd  = (( ref_drydp['DryDep_PAN'] + ref_drydp['DryDep_BZPAN'] + ref_drydp['DryDep_MPAN'] )
                                         * hemco['AREA'] * 1.0e4 * 365. * 24. * 3600 / avo * 14.007 * 1e-12).sum().values
    N2O5_dd = ( ref_drydp['DryDep_N2O5'] * hemco['AREA'] * 1.0e4 * 365. * 24. * 3600 / avo * 14.007 * 1e-9).sum().values 
    HNO2_dd = ( ref_drydp['DryDep_HNO2'] * hemco['AREA'] * 1.0e4 * 365. * 24. * 3600 / avo * 14.007 * 1e-9).sum().values
    NIT_dd  = (( ref_drydp['DryDep_NIT'] + ref_drydp['DryDep_NITs'] + ref_drydp['DryDep_NITD1'] +
                 ref_drydp['DryDep_NITD2'] + ref_drydp['DryDep_NITD3'] + ref_drydp['DryDep_NITD4'] )
                                         * hemco['AREA'] * 1.0e4 * 365. * 24. * 3600 / avo * 14.007 * 1e-12).sum().values
    Hal_dd  = (( ref_drydp['DryDep_ClNO3'] + ref_drydp['DryDep_BrNO3'] )
                                         * hemco['AREA'] * 1.0e4 * 365. * 24. * 3600 / avo * 14.007 * 1e-9).sum().values



if __name__=="__main__":
    main()
