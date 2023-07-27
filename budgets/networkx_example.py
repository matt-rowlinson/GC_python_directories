import networkx as nx
import matplotlib.pyplot as plt
import xarray as xr

pos = {
       'pNO$_3^-$'  : (-10,0),
       'NO'         : (26,-28),
       'NO$_2$'     : (-8,-29),
       'N$_2$O$_5$' : (2,29),
       'HNO$_3$'    : (28,22),
       'PAN'        : (-28, 28),
       'HONO'       : (12,2),
       'NO$_3$'        : (-28,-25)}

options = {
            "node_size": 7e3,
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
          'NO$_3$'    :'NO$_3$'}

gas_phase_edge_list = { ('pNO$_3^-$','HONO'),
                        ('pNO$_3^-$','N$_2$O$_5$'),
                        ('pNO$_3^-$','NO$_2$'),
                        ('HONO','NO'),
                        ('HNO$_3$','HONO'),
                        ('NO', 'HONO'),
                        ('NO', 'NO$_2$'),
                        ('NO$_2$','NO$_3$'),
                        ('NO$_2$','PAN'),
                        ('NO$_2$','HONO'),
                        ('NO$_2$','NO')}

het_edge_list = { ('HNO$_3$','NO$_2$'),
                  ('HNO$_3$','pNO$_3^-$')}

def draw_my_graph(pos,
                  options,
                  Labels,
                  gas_phase_edge_list,
                  het_edge_list,
                  edge_labels
                  ):
    f=plt.figure(figsize=(12,8))
    G = nx.Graph()

    nx.draw_networkx_nodes( G, pos, nodelist=Labels.keys(), **options)
    nx.draw_networkx_labels( G, pos, Labels, font_size=26 )
    nx.draw_networkx_edges( G,pos,edgelist=gas_phase_edge_list,
                                  arrows=True, arrowstyle='-|>',width=5, node_size=7e3,                              
                                  connectionstyle='arc3,rad=-0.1') 
    nx.draw_networkx_edges( G,pos,edgelist=het_edge_list,edge_color='r',
                                  arrows=True, arrowstyle='->',style='dashed',width=2.5, node_size=7e3,  
                                  connectionstyle='arc3,rad=0.31')

    nx.draw_networkx_edge_labels( G,pos,edge_labels,
                                  verticalalignment='center', font_size=16
                                  )

    # Set margins for the axes so that nodes aren't clipped
    ax = plt.gca()
    ax.set_xlim( -35,35 )
    ax.set_ylim( -35,35 )
    #ax.margins(0.01)
    plt.tight_layout()
    plt.axis("off")
    plt.savefig('plots/my-budget-plot.png')
    return 

####------------------------
def main():
    ## Read GEOS-Chem tagged run output
    ds0 = 


    a=0.0
    b=2.7
    edge_labels = { ('HONO','NO'):f"{b}\n({a})",
                    ('NO$_2$','PAN'):102.2     
                    }
    draw_my_graph(pos,options,Labels,gas_phase_edge_list,het_edge_list,edge_labels)

if __name__=="__main__":
    main()
