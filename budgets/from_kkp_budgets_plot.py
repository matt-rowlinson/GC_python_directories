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
       'PAN'                : (-28, 28),
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
    tlist = [ "SpeciesConcVV_" + s for s in tracers]
    
    ds = ref_spc[tlist].to_array().sum("variable")
    ds = ds * ref_met['Met_AD'].values * 1.0e3 * MW / MW_air * 1.0e-9
    Trop_Mask = (ref_met['Met_PMID'] > ref_met['Met_TropP'])#.mean(dim='time')
    ds0 = ds.where(Trop_Mask)

    ds = dev_spc[tlist].to_array().sum("variable")
    ds = ds * ref_met['Met_AD'].values * 1.0e3 * MW / MW_air * 1.0e-9
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
    ref_met = get_data_as_xr('new_base_4x5_tracers',dates=["StateMet.2019"]).mean(dim='time')
    ref_spc = get_data_as_xr('new_base_4x5_tracers',dates=["SpeciesConc.2019"]).mean(dim='time')

    dev_met = get_data_as_xr('langmuir_4x5_tracers',dates=["StateMet.2019"])
    dev_spc = get_data_as_xr('langmuir_4x5_tracers',dates=["SpeciesConc.2019"])


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
        label_rxns( ref_met, ref_spc,
                    dev_met, dev_spc,
                    MW,
                    G, edge, tracers, font_color='k' )
    #----------------------------
    ax = plt.gca()
    ax.set_xlim( -35,35 )
    ax.set_ylim( -35,35 )
    plt.tight_layout()
    plt.axis("off")
    plt.savefig('plots/KPP-filled-budget-plot.png')
    plt.close()
    

if __name__=="__main__":
    main()
