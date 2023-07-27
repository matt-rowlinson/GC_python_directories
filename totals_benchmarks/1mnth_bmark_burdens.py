import pandas as pd
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("agg")
from matplotlib.pyplot import figure
figure(figsize=(8, 8), dpi=80)

path='./*.0'

versions=[] ; over=[]
for version in sorted(glob.glob(path)):
    versions.append( version.replace('./','' ) )
    infile=f'{version}/GlobalMass_Trop.txt'
    print( infile )
    df = pd.read_csv( infile, header=4,sep='\s{2,}')
    over.append( df.Dev )
    
df = pd.concat( over, axis=1 )
df.columns= versions 
Is=[]
for i, s in df.iterrows():
    plt.plot( s )
    i = str(i).replace('(','').replace(')','').replace(':','')\
            .replace(',','').replace("'",'').replace(' ','')
    Is.append( i )
    plt.ylabel( f'{i} (Gg)' )
    plt.savefig( f'1month_plots/{i}_burdens.png' )
    plt.close()

df.index = Is
NOx = df.loc['NO2'] + df.loc['NO']
plt.plot( NOx )
plt.ylabel( 'NOx (Gg)' )
plt.savefig( f'1month_plots/NOx_burdens.png' )
plt.close()

# NOy, Ox, Bry, Cly, Iy
NOy = df.loc['NO'] + df.loc['NO2'] + df.loc['N2O5'] + df.loc['ClNO2'] + df.loc['HNO3'] + df.loc['NIT'] + df.loc['NITs'] + \
      df.loc['HNO2'] + df.loc['NH3'] + df.loc['PAN'] + df.loc['HNO4'] 
plt.plot( NOy )
plt.ylabel( 'NOy (Gg)' )
plt.savefig( f'1month_plots/NOy_burdens.png' )
plt.close()

Bry = df.loc['Br'] + (df.loc['Br2']*2) + df.loc['HOBr'] + df.loc['BrO'] + df.loc['HBr'] + df.loc['BrNO2'] + df.loc['BrNO3'] + \
      df.loc['IBr'] + df.loc['BrCl']
plt.plot( Bry )
plt.ylabel( 'Bry (Gg)' )
plt.savefig( f'1month_plots/Bry_burdens.png' )
plt.close()

Cly = df.loc['Cl'] + (df.loc['Cl2']*2) + df.loc['HOCl'] + df.loc['ClO'] + df.loc['HCl'] + df.loc['ClNO2'] + df.loc['ClNO3'] + \
      df.loc['ICl'] + df.loc['BrCl'] + df.loc['ClOO'] + df.loc['OClO'] + df.loc['Cl2O2']
plt.plot( Cly )
plt.ylabel( 'Cly (Gg)' )
plt.savefig( f'1month_plots/Cly_burdens.png' )
plt.close()

Iy = df.loc['I'] + (df.loc['I2']*2) + df.loc['HOI'] + df.loc['IO'] + df.loc['OIO'] + df.loc['HI']
plt.plot( Iy )
plt.ylabel( 'Iy (Gg)' )
plt.savefig( f'1month_plots/Iy_burdens.png' )
plt.close()

Ox = df.loc['BrNO2']+ (df.loc['BrNO3']*2)+ df.loc['BrO']+ (df.loc['Cl2O2']*2)+ df.loc['ClNO2']+ (df.loc['ClNO3']*2)+ df.loc['ClO']\
        + df.loc['ETHLN']+ df.loc['ETNO3']+ df.loc['HNO3']+ df.loc['HNO4']+ df.loc['HOBr']+ df.loc['HOCl']+ df.loc['HOI']\
        + (df.loc['I2O2']*2)+ (df.loc['I2O3']*3)+ (df.loc['I2O4']*4)+ df.loc['ICN']+ df.loc['ICNOO']+ df.loc['IDHNBOO']\
        + df.loc['IDHNDOO1']+ df.loc['IDHNDOO2']+ (df.loc['IDN']*2)+ (df.loc['IDNOO']*2)+ df.loc['IHN1']+ df.loc['IHN2']\
        + df.loc['IHN3']+ df.loc['IHN4']+ df.loc['IHPNBOO']+ df.loc['IHPNDOO']+ df.loc['INA']+ df.loc['INO2B']+ df.loc['INO2D']\
        + df.loc['INPB']+ df.loc['INPD']+ df.loc['IO']+ df.loc['IONO']+ (df.loc['IONO2']*2)+ df.loc['ISOPNOO1']+ df.loc['ISOPNOO2']\
        + df.loc['ITCN']+ df.loc['ITHN']+ df.loc['MACRNO2']+ df.loc['MCRHN']+ df.loc['MCRHNB']+ df.loc['MENO3']+ df.loc['MPAN']\
        + df.loc['MPN']+ df.loc['MVKN']+ (df.loc['N2O5']*3)+ df.loc['NO2']+ (df.loc['NO3']*2)+ df.loc['NPRNO3']+ df.loc['O']\
        + df.loc['O1D']+ df.loc['O3']+ (df.loc['OClO']*2)+ (df.loc['OIO']*2)+ df.loc['OLND']+ df.loc['OLNN']+ df.loc['PAN']\
        +df.loc['PPN']+ df.loc['PRN1']+ df.loc['PROPNN']+ df.loc['PRPN']+ df.loc['R4N1']+ df.loc['R4N2']
plt.plot( Ox )
plt.ylabel( 'Ox (Gg)' )
plt.savefig( f'1month_plots/Ox_burdens.png' )
plt.close()


plt.plot( df.loc['O3'], label='O3' )
plt.plot( Ox, label='Ox' )
plt.plot( Bry + Cly + Iy, label='Halogens' )
plt.plot( NOy, label='NOy' )
plt.legend()
plt.ylabel( 'Gg' )
plt.savefig( f'1month_plots/together_burdens.png' )
plt.close()
