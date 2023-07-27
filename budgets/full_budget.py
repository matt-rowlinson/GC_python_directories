import networkx as nx
import matplotlib.pyplot as plt
import xarray as xr
import os
import sys
import numpy as np
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

edge_list = {           ('NO', 'NO2')               :['NG001','NG003','NG004','NG007'],
                        ('NO2','NO')                :['NG002','NG011','OP001'],
                        ('HONO', 'NO2')             :['NG005'],
                        ('NO3', 'NO2')              :['NG006','NG007','NG008','OP004'],
                        ('NO2', 'N2O5')             :['NG009'],
                        ('NO3', 'N2O5')             :['NG009'],
                        ('N2O5', 'NO2')             :['NG010','OP006'],
                        ('N2O5', 'NO3')             :['NG010','OP006'],
                        ('NO3', 'NO')               :['NG011','OP005'],
                        ('NO2','PAN')               :['NG012'],
                        ('PAN','NO2')               :['NG013','OP007'],
                        ('NO', 'HONO')              :['NG014'],
                        ('HNO3', 'NO3')             :['NG015'],
                        ('pNO3','HONO')             :['NP001','NP002','NP003'],
                        ('pNO3','NO2')              :['NP004','NP005','NP006'],
                        ('HNO3','NO2')              :['OP002'],
                        ('HONO','NO')               :['OP003'],
                        ('PAN','NO3')               :['OP007'],
                        ('NO2','HONO')              :['HH001'],
                        ('NO2','HNO3')              :['HH001'],
                        ('NO3','HNO3')              :['HH002'],
                        ('N2O5','HNO3')             :['HH003'],
                        ('Halogen nitrates','HNO3') :['HH004','HH005','HH006'],
                        ('Halogen nitrates','HONO') :['HN001','HN002']
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

def draw_my_graph(pos,
                  options,
                  Labels,
                  edge_list,
                  rad1=-0.1, rad2=0.3
                  ):
    f=plt.figure(figsize=(12,8))
    G = nx.Graph()

    nx.draw_networkx_nodes( G, pos, nodelist=Labels.keys(), **options)
    nx.draw_networkx_labels( G, pos, Labels, font_size=26 )
    nx.draw_networkx_edges( G,pos,edgelist=edge_list.keys(),
                                  arrows=True, arrowstyle='-|>',width=5, node_size=7e3,                              
                                  connectionstyle=f'arc3,rad={rad1}') 
    return G 

def label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                MW,
                G, edge, tracers,label_pos=.7,font_color='k',rotate=False):
    avo=6.022e23
    MW_air = 28.96
    tlist = [ "Prod_" + s for s in tracers]
    
    #ds = ref_spc[tlist].to_array().sum("variable")
    #ds = ds * ref_met['Met_AD'].values * 1.0e3 * MW / MW_air * 1.0e-9
    #Trop_Mask = (ref_met['Met_PMID'] > ref_met['Met_TropP'])#.mean(dim='time')
    #ds0 = ds.where(Trop_Mask)

    ds = dev_spc[tlist].to_array().sum("variable")
    #ds = ds * ref_met['Met_AD'].values * 1.0e3 * MW / MW_air * 1.0e-9
    ds = ds * dev_met['Met_AIRVOL'] *  1.0e6   * 365.0 * 24.0 * 3600.0 / avo * MW * 1.0e-12

    Trop_Mask = (dev_met['Met_PMID'] > dev_met['Met_TropP'])#.mean(dim='time')
    ds1 = ds.where(Trop_Mask)

    label = {edge: '%.2f' %ds1.sum().values} ## Single value only
    #label = {edge: '%.2f\n$^{(%.2f)}$' %(ds1.sum().values,ds0.sum().values)} ## Base value in brackets 

    nx.draw_networkx_edge_labels( G,pos,label,horizontalalignment='center',
            verticalalignment='center', font_size=12,label_pos=label_pos,font_color=font_color, rotate=rotate)
    return

####------------------------
def main():
    
    ## Read GEOS-Chem tagged run output
    ref_met = get_data_as_xr('langmuir_4x5_tracers',dates=["StateMet.2019"]).mean(dim='time')
    ref_spc = get_data_as_xr('langmuir_4x5_tracers',dates=["SpeciesConc.2019"]).mean(dim='time')
    ref_pro = get_data_as_xr('langmuir_4x5_tracers',dates=["ProdLoss.20190810_0200"]).mean(dim='time')

    dev_met = get_data_as_xr('langmuir_4x5_tracers',dates=["StateMet.2019"])
    dev_spc = get_data_as_xr('langmuir_4x5_tracers',dates=["SpeciesConc.2019"])
    dev_pro = get_data_as_xr('langmuir_4x5_tracers',dates=["ProdLoss.20190810_0200"])

    # Draw the networkx graph
    G = draw_my_graph(pos,options,Labels,edge_list,
                        rad1=0., rad2=0.)

    # Label each pathway one by one using the dictionary of pathways and tracers
    for edge, tracers in edge_list.items():
        print( edge, tracers )
        if edge[1]=="N2O5":
            MW=14.007*2
        else:
            MW=14.007
        label_rxns( ref_met, ref_pro,
                    dev_met, dev_pro,
                    MW,
                    G, edge, tracers, font_color='k' )
    
    avo=6.022e23
    hemco = get_data_as_xr('langmuir_RH_904x5_tracers',dates=["HEMCO_diag"]).mean(dim='time')
    ref_prodl = get_data_as_xr('langmuir_RH_904x5_tracers',dates=["ProdLoss.2019"]).mean(dim='time')
    ref_drydp = get_data_as_xr('langmuir_RH_904x5_tracers',dates=["DryDep.2019"]).mean(dim='time')

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


    #----------------------------
    ax = plt.gca()
    
    ## Add emissions to plot
    NOx     = ax.text(10. , -45, "NO$_x$",fontsize=15,weight='bold',ha="left",)
    NOx_em  = ax.text(11, -37, "%.f TgN" %NOx_em, ha="left", va="center", rotation=90, size=15,
              bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightblue", ec="steelblue", lw=2))
    HNO3    = ax.text(16, -45, "HNO$_3$",fontsize=15,weight='bold' ,ha="left",)
    HNO3_em = ax.text(17, -37, "%.f TgN" %HNO3_em, ha="left", va="center", rotation=90, size=15,
              bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightblue", ec="steelblue", lw=2))

    ## Add dry deposition to plot
    NO2_dd  = ax.text(-10, -38, f"%.f TgN" %NO2_dd, ha="center", va="center", rotation=45, size=15,
              bbox=dict(boxstyle="larrow,pad=0.3",fc="lightyellow", ec="yellow", lw=2))
    HNO2_dd = ax.text(4, 20, f"%.f GgN" %HNO2_dd, ha="center", va="center", rotation=90, size=15,
               bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightyellow", ec="yellow", lw=2))
    HNO3_dd = ax.text(30, 10, f"%.f TgN" %HNO3_dd, ha="center", va="center", rotation=315, size=15,
              bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightyellow", ec="yellow", lw=2))
    PAN_dd  = ax.text(-28, 33, f"%.f TgN" %PAN_dd, ha="center", va="center", rotation=90, size=15,
              bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightyellow", ec="yellow", lw=2))
    NIT_dd  = ax.text(-7, 8, f"%.f TgN" %NIT_dd, ha="center", va="center", rotation=90, size=15,
              bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightyellow", ec="yellow", lw=2))
    N2O5_dd = ax.text(-16, 32, f"%.f GgN" %N2O5_dd, ha="center", va="center", rotation=90, size=15,
              bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightyellow", ec="yellow", lw=2))
    Hal_dd  = ax.text(26, 33, f"%.f GgN" %Hal_dd, ha="center", va="center", rotation=25, size=15,
              bbox=dict(boxstyle="rarrow,pad=0.3",fc="lightyellow", ec="yellow", lw=2))


    ax.set_xlim( -40,40 )
    ax.set_ylim( -45,40 )
    plt.tight_layout()
    plt.axis("off")
    plt.savefig('plots/langmuir_RH_90full-budget-plot.png')
    plt.close()

if __name__=="__main__":
    main()
