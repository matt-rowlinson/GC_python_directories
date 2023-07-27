#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
plt.style.use('seaborn-darkgrid')

def main():
    df0 = pd.read_csv('csvs/emep-ceds_ratios_2017.csv', index_col=0)[['NAEI','v2021-06']]
    df1 = pd.read_csv('csvs/nei-ceds_ratios_2017.csv', index_col=0)[['NEI','v2021-06']]
    df2 = pd.read_csv('csvs/meic-ceds_ratios_2017.csv', index_col=0)[['MEIC','v2021-06']]
    index = ['HCHO','Ethanol','Ethane','Propane','Propene','Ketones','ALK4','Benzene','Toluene','Xylene']
    print( df1 )
    #df1 = df1[['NEI','v2021-06']]
    #print( df1 )

    c = ['#66c2a5','#fc8d62','#8da0cb']
    f,(ax1,ax2,ax3)=plt.subplots(3,1, figsize=(8,8))
    df0.plot.bar(rot=0, ax=ax1, width=.4, color=c)
    df1.plot.bar(rot=0, ax=ax2, width=.4, color=c)
    df2.plot.bar(rot=0, ax=ax3, width=.4, color=c)

    ax1.set_ylabel('VOC : Total NMVOCs')
    ax2.set_ylabel('VOC : Total NMVOCs')
    ax3.set_ylabel('VOC : Total NMVOCs')

    ax1.set_xticks(range(0,10), index)
    ax2.set_xticks(range(0,10), index)
    ax3.set_xticks(range(0,10), index)


    plt.tight_layout()
    plt.savefig(f'plots/figure_1.png' )
    plt.close()

if __name__ == "__main__":
    main()
