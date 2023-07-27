import networkx as nx
import matplotlib.pyplot as plt
import xarray as xr
import os
import sys
import numpy as np
#--------------------------------
pos = {
       'pNO$_3^-$'  : (-7,0),
       'NO'         : (26,-20),
       'NO$_2$'     : (-8,-31),
       'N$_2$O$_5$' : (-16,23),
       'HNO$_3$'    : (28,17),
       'PAN'        : (-28, 28),
       'HONO'       : (12,7),
       'Halogen\nnitrates'       : (16,30),
       'NO$_3$'        : (-28,-25)}

options = {
            "node_size": 7e2,
            "node_color": "white",
            "edgecolors": "white",
            "alpha":0.9,
            "node_shape":'o'}
        
Labels = {'NO'        :'NO',
          'HONO'      :'HONO',
          'PAN'       :'PAN',
          'pNO$_3^-$' :'pNO$_3^-$',
          'N$_2$O$_5$':'N$_2$O$_5$',
          'NO$_2$'    :'NO$_2$',
          'HNO$_3$'   :'HNO$_3$',
          'Halogen\nnitrates'    :'(Br,Cl,I)NO$_3$',
          'NO$_3$'    :'NO$_3$'}

gas_phase_edge_list = { ('NO', 'NO$_2$')        :['NG001','NG003','NG004','NG007'],
                        ('NO$_2$','NO')         :['NG002','NG011'],
                        ('HONO', 'NO$_2$')      :['NG005'],
                        ('NO$_3$', 'NO$_2$')    :['NG006','NG007','NG008'],
                        ('NO$_2$', 'N$_2$O$_5$'):['NG009'],
                        ('NO$_3$', 'N$_2$O$_5$'):['NG009'],
                        ('N$_2$O$_5$', 'NO$_2$'):['NG010'],
                        ('N$_2$O$_5$', 'NO$_3$'):['NG010'],
                        ('NO$_3$', 'NO')        :['NG011'],
                        ('NO$_2$','PAN')        :['NG012'],
                        ('PAN','NO$_2$')        :['NG013'],
                        ('NO', 'HONO')          :['NG014'],
                        ('HNO$_3$', 'NO$_3$')   :['NG015'],
                        }

photolysis_edge_list = {('pNO$_3^-$','HONO')    :['NP001','NP002','NP003'],
                        ('pNO$_3^-$','NO$_2$')      :['NP004','NP005','NP006'],
                        ('NO$_2$','NO')             :['OP001'],
                        ('HNO$_3$','NO$_2$')        :['OP002'],
                        ('HONO','NO')               :['OP003'],
                        ('NO$_3$','NO$_2$')         :['OP004'],
                        ('NO$_3$','NO')             :['OP005'],
                        ('N$_2$O$_5$','NO$_2$')     :['OP006'],
                        ('N$_2$O$_5$','NO$_3$')     :['OP006'],
                        ('PAN','NO$_2$')            :['OP007'],
                        ('PAN','NO$_3$')            :['OP007']
                        }          

het_edge_list = {('NO$_2$','HONO')                :['HH001'],
                 ('NO$_2$','HNO$_3$')             :['HH001'],
                 ('NO$_3$','HNO$_3$')             :['HH002'],
                 ('N$_2$O$_5$','HNO$_3$')         :['HH003'],
                 ('Halogen\nnitrates','HNO$_3$')  :['HH004','HH005','HH006'],
                 ('Halogen\nnitrates','HONO')     :['HN001','HN002']
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
                  gas_phase_edge_list,
                  photolysis_edge_list,
                  het_edge_list,
                  rad1=-0.1, rad2=0.3
                  ):
    f=plt.figure(figsize=(12,8))
    G = nx.Graph()

    nx.draw_networkx_nodes( G, pos, nodelist=Labels.keys(), **options)
    nx.draw_networkx_labels( G, pos, Labels, font_size=26 )
    nx.draw_networkx_edges( G,pos,edgelist=gas_phase_edge_list.keys(),
                                  arrows=True, arrowstyle='-|>',width=5, node_size=7e3,                              
                                  connectionstyle=f'arc3,rad={rad1}') 
    nx.draw_networkx_edges( G,pos,edgelist=photolysis_edge_list,edge_color='y',
                                  arrows=True, arrowstyle='->',width=2.5, node_size=7e3,  
                                  connectionstyle=f'arc3,rad={rad1}')
    nx.draw_networkx_edges( G,pos,edgelist=het_edge_list,edge_color='r',
                                  arrows=True, arrowstyle='->',style='dashed',width=2.5, node_size=7e3,  
                                  connectionstyle=f'arc3,rad={rad2}')

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

    #label = {edge: '%.2f' %ds1.sum().values} ## Single value only
    label = {edge: '%.2f\n$^{(%.2f)}$' %(ds1.sum().values,ds0.sum().values)} ## Base value in brackets 

    nx.draw_networkx_edge_labels( G,pos,label,horizontalalignment='center',
            verticalalignment='center', font_size=12,label_pos=label_pos,font_color=font_color, rotate=rotate)
    return

####------------------------
def main():
    avo=6.022e23

    ## Read GEOS-Chem tagged run output
    ref_met = get_data_as_xr('new_base_4x5_tracers',dates=["StateMet.2019"]).mean(dim='time')
    ref_spc = get_data_as_xr('new_base_4x5_tracers',dates=["SpeciesConc.2019"]).mean(dim='time')

    dev_met = get_data_as_xr('langmuir_4x5_tracers',dates=["StateMet.2019"])
    dev_spc = get_data_as_xr('langmuir_4x5_tracers',dates=["SpeciesConc.2019"])


    # Draw the networkx graph
    G = draw_my_graph(pos,options,Labels,gas_phase_edge_list,photolysis_edge_list,het_edge_list,
                        rad1=0., rad2=0.)

    # Label each pathway one by one using the dictionary of pathways and tracers
    for edge, tracers in gas_phase_edge_list.items():
        print( edge, tracers )
        label_rxns( ref_met, ref_spc,
                    dev_met, dev_spc,
                    14.007,
                    G, edge, tracers, font_color='k' )

    for edge, tracers in photolysis_edge_list.items():
        print( edge, tracers )
        label_rxns( ref_met, ref_spc,
                    dev_met, dev_spc,
                    14.007,
                    G, edge, tracers, font_color='k' )

    for edge, tracers in het_edge_list.items():
        print( edge, tracers )
        label_rxns( ref_met, ref_spc,
                    dev_met, dev_spc,
                    14.007,
                    G, edge, tracers, font_color='k' )
    #sys.exit()

    #----------------------------
    ax = plt.gca()
    ax.set_xlim( -35,35 )
    ax.set_ylim( -35,35 )
    plt.tight_layout()
    plt.axis("off")
    plt.savefig('plots/color-filled-budget-plot.png')
    plt.close()
    

if __name__=="__main__":
    main()
