import networkx as nx
import matplotlib.pyplot as plt
import xarray as xr
import os
import sys
import numpy as np
#--------------------------------
pos = {
       'pNO$_3^-$'  : (-10,0),
       'NO'         : (26,-28),
       'NO$_2$'     : (-8,-29),
       'N$_2$O$_5$' : (-12,29),
       'HNO$_3$'    : (28,22),
       'PAN'        : (-28, 28),
       'HONO'       : (12,2),
       'Halogen\nnitrates'       : (6,30),
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

gas_phase_edge_list = { ('pNO$_3^-$','HONO'),
                        ('pNO$_3^-$','NO$_2$'),
                        ('NO$_3$','N$_2$O$_5$'),
                        ('HONO','NO'),
                        ('HNO$_3$','HONO'),
                        ('NO', 'HONO'),
                        ('NO', 'NO$_2$'),
                        ('NO$_2$','NO$_3$'),
                        ('NO$_3$','NO$_2$'),
                        ('NO$_2$','PAN'),
                        ('HONO','NO$_2$'),
                        ('NO$_2$','NO')}

het_edge_list = { ('HNO$_3$','NO$_2$'),
                  ('HNO$_3$','pNO$_3^-$'),
                  ('Halogen\nnitrates','HNO$_3$'),
                  ('Halogen\nnitrates','HONO'),
                  ('NO$_2$','HONO')}
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
                  het_edge_list,
                  rad1=-0.1, rad2=0.3
                  ):
    f=plt.figure(figsize=(12,8))
    G = nx.Graph()

    nx.draw_networkx_nodes( G, pos, nodelist=Labels.keys(), **options)
    nx.draw_networkx_labels( G, pos, Labels, font_size=26 )
    nx.draw_networkx_edges( G,pos,edgelist=gas_phase_edge_list,
                                  arrows=True, arrowstyle='-|>',width=5, node_size=7e3,                              
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
    G = draw_my_graph(pos,options,Labels,gas_phase_edge_list,het_edge_list,
                        rad1=0., rad2=0.3)

    # Label each pathway one by one using a list of tracers 
    edge    = ('pNO$_3^-$','HONO')
    tracers = ['NP001','NP002','NP003']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )

    edge    = ('pNO$_3^-$','NO$_2$')
    tracers = ['NP004','NP005','NP006']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )

    edge    = ('NO$_2$','NO')
    tracers = ['OP001','NG002']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )

    edge    = ('NO', 'NO$_2$')
    tracers = ['NG001','NG003','NG007']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )  
    
    edge    = ('NO$_3$','N$_2$O$_5$')
    tracers = ['NG009']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007*2.,
                G, edge, tracers )

    edge    = ('NO$_2$','PAN')
    tracers = ['NG012']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )

    edge    = ('PAN','NO$_2$')
    tracers = ['NG013']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )
    
    edge    = ('HONO','NO')
    tracers = ['OP003']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )
    
    edge    = ('HONO','NO$_2$')
    tracers = ['NG005']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )

    edge    = ('NO','HONO')
    tracers = ['NG014']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )

    edge    = ('NO$_2$','HONO')
    tracers = ['HH001']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers, font_color='r' )

    edge    = ('Halogen\nnitrates','HNO$_3$')
    tracers = ['HH004','HH005','HH006']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers, font_color='r' )

    edge    = ('Halogen\nnitrates','HONO')
    tracers = ['HN001','HN002']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers, font_color='r' )

    edge    = ('N$_2$O$_5$','HNO$_3$')
    tracers = ['HH003']
    label_rxns( ref_met, ref_spc,
                dev_met, dev_spc,
                14.007,
                G, edge, tracers )

    #----------------------------
    ax = plt.gca()
    ax.set_xlim( -35,35 )
    ax.set_ylim( -35,35 )
    plt.tight_layout()
    plt.axis("off")
    plt.savefig('plots/filled-budget-plot.png')
    plt.close()
    

if __name__=="__main__":
    main()
