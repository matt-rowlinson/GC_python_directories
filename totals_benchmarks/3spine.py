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
    hold=[]
    for infile in sorted(glob.glob(f'{version}/*')):
        print( infile )
        df = pd.read_csv( infile, header=4,sep='\s{2,}')
        hold.append( df.Dev )
    
    df = pd.concat( hold, axis=1 )
    df = df.mean(axis=1)
    over.append( df )

df = pd.concat( over, axis=1 )
df.columns= versions 
Is=[]
for i, s in df.iterrows():
    i = str(i).replace('(','').replace(')','').replace(':','')\
            .replace(',','').replace("'",'').replace(' ','')
    Is.append( i )

df.index = Is
NOx = df.loc['NO2'] + df.loc['NO']
NOy = df.loc['NO'] + df.loc['NO2'] + df.loc['N2O5'] + df.loc['ClNO2'] + df.loc['HNO3'] + df.loc['NIT'] + df.loc['NITs'] + \
      df.loc['HNO2'] + df.loc['NH3'] + df.loc['PAN'] + df.loc['HNO4'] 
Bry = df.loc['Br'] + (df.loc['Br2']*2) + df.loc['HOBr'] + df.loc['BrO'] + df.loc['HBr'] + df.loc['BrNO2'] + df.loc['BrNO3'] + \
      df.loc['IBr'] + df.loc['BrCl']
Cly = df.loc['Cl'] + (df.loc['Cl2']*2) + df.loc['HOCl'] + df.loc['ClO'] + df.loc['HCl'] + df.loc['ClNO2'] + df.loc['ClNO3'] + \
      df.loc['ICl'] + df.loc['BrCl'] + df.loc['ClOO'] + df.loc['OClO'] + df.loc['Cl2O2']
Iy = df.loc['I'] + (df.loc['I2']*2) + df.loc['HOI'] + df.loc['IO'] + df.loc['OIO'] + df.loc['HI']
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

'''
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

host = host_subplot(111, axes_class=AA.Axes)
f, ax = plt.subplots(1)
plt.subplots_adjust(right=0.75)

par1 = host.twinx()
par2 = host.twinx()
par1 = ax.twinx()
par2 = ax.twinx()

offset = 60
#new_fixed_axis = par2.get_grid_helper().new_fixed_axis
#par2.axis["right"] = new_fixed_axis(loc="right",
#        axes=par2,
#        offset=(offset, 0))

#par2.axis["right"].toggle(all=True)

#host.set_xlim(0, 2)
#host.set_ylim(0, 2)

host.set_xlabel("Version")
host.set_ylabel("Ox")
par1.set_ylabel("NOy")
par2.set_ylabel("Halogens")

#p1, = host.plot( df.columns, df.loc['O3'], label="Ox")
p1, = ax.plot( df.columns, df.loc['O3'], label="Ox")


print( df.loc['O3'])
print( NOy )

p2, = par1.plot( df.columns, NOy, label="NOy")
p3, = par2.plot( df.columns, Bry + Cly + Iy, label="Halogens")

#par1.set_ylim(0, 4)
#par2.set_ylim(1, 65)

host.legend()

#host.axis["left"].label.set_color(p1.get_color())
#par1.axis["right"].label.set_color(p2.get_color())
#par2.axis["right"].label.set_color(p3.get_color())

plt.draw()
plt.savefig('plots/together_burdens.png')
plt.close()
'''

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)

twin1 = ax.twinx()
twin2 = ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
twin2.spines['right'].set_position(("axes", 1.2))
p1, = ax.plot(df.loc['O3'], "b", label="O3")
p2, = twin1.plot(NOx, "r", label="NOx")
p3, = twin2.plot(Bry + Cly + Iy, "g", label="Halogens")

ax.set_xlabel("Version")
ax.set_ylabel("O3 (Gg)")
twin1.set_ylabel("NOx (Gg)")
twin2.set_ylabel("Halogens (Gg)")

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())

tkw = dict(size=4, width=1.5)
ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)

ax.legend(ncol=3, handles=[p1, p2, p3])

plt.savefig('plots/together_burdens_NOx.png')
plt.close()



fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)

twin1 = ax.twinx()
twin2 = ax.twinx()

# Offset the right spine of twin2.  The ticks and label have already been
# placed on the right by twinx above.
twin2.spines['right'].set_position(("axes", 1.2))
p1, = ax.plot( Bry, "b", label="Bry")
p2, = twin1.plot(df.loc['SALCAL'], "r", label="SALCAL")
p3, = twin2.plot(df.loc['SALAAL'], "g", label="SALAAL")

ax.set_xlabel("Version")
ax.set_ylabel("Bry (Gg)")
twin1.set_ylabel("SALCAL (Gg)")
twin2.set_ylabel("SALAAL (Gg)")

ax.yaxis.label.set_color(p1.get_color())
twin1.yaxis.label.set_color(p2.get_color())
twin2.yaxis.label.set_color(p3.get_color())

tkw = dict(size=4, width=1.5)
ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)
ax.tick_params(axis='x', **tkw)

ax.legend(ncol=3, handles=[p1, p2, p3])

plt.savefig('plots/alkilinity_burdens_Bry.png')
plt.close()
