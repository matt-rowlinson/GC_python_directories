


import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use("agg")
plt.style.use('seaborn-darkgrid')

f=plt.figure(figsize=(10,5))
ax1 = f.add_subplot(1,2,1)
# Ethane emissions
xticks=['Global','USA','Europe','China','NH','SH']
ceds  = [  7.06,  0.33, 0.71, 1.2, 5.67, 1.39]
ceds21= [ 4.94, 0.25, 0.25+0.41, 2.47, 4.12, 0.82]
#edgar = [ 4.855, 0.29, 0.50, 0.65 ]
me    = [ 10.96, 2.03, 1.55, 2.88, 9.57, 1.39 ]
me2   = [ 8.97 , 1.78, 0.51+0.82, 4.45, 8.14, 0.82 ]
#dalso = [ 13.60, np.nan, np.nan, np.nan ]

cs = ['#1b9e77','#d95f02','#7570b3','#e7298a']
ax1.plot( [xticks[0],xticks[0]], [3.36, 10.],marker='_', solid_capstyle='butt', lw=6,c='y', alpha=.3, label='Bottom up estimates')
ax1.scatter( xticks, ceds , s=42, c=cs[0], marker='o', label='CEDS 2020-08' )
#ax1.scatter( xticks, edgar, s=42, c=cs[1], marker='^', label='EGDAR')
#ax1.scatter( xticks, dalso, s=52, c=cs[3], marker='+', label='Dalsøren et al. (2019)')
ax1.scatter( xticks, me   , s=82, c=cs[2], marker='*', label='This study')
ax1.scatter( xticks, ceds21, s=102,c='r',   marker='o', label='CEDS 2021-06' )
ax1.scatter( xticks, me2   , s=102,c='gold',marker='*', label='NEW' )
ax1.set_ylabel('Tg yr$^{-1}$')
ax1.set_title('Anthropogenic ethane emissions')
plt.legend()

ax2 = f.add_subplot(1,2,2)
## Propane emissions
ceds  = [ 6.40, 0.36, 0.83, 1.2, 5.03, 1.37 ]
ceds21= [ 4.13, 0.33, 0.27+0.44, 1.66, 3.51, 0.62 ]
#edgar = [ 6.22, 0.2379, 0.411, 1.4 ]
me    = [ 8.10, 1.21, 1.39, 1.88, 6.73, 1.37 ]
me2   = [ 6.49, 1.71, 0.54+0.88, 2.19, 5.87, 0.62 ]
#dalso = [ 16.10, np.nan, np.nan, np.nan ]

ax2.plot( [xticks[0],xticks[0]], [3., 10.],marker='_', solid_capstyle='butt', lw=6,c='y', alpha=.3, label='Bottom up estimates')
ax2.scatter( xticks, ceds , s=42, c=cs[0], marker='o', label='CEDS 2020-08' )
#ax2.scatter( xticks, edgar, s=42, c=cs[1], marker='^', label='EDGAR')
#ax2.scatter( xticks, dalso, s=52, c=cs[3], marker='+', label='Dalsøren et al. (2019)')
ax2.scatter( xticks, me   , s=82, c=cs[2], marker='*', label='This study')
ax2.scatter( xticks, ceds21,s=102,c='r',   marker='o', label='CEDS 2021-06' )
ax2.scatter( xticks, me2  , s=102,c='gold',marker='*', label='NEW')
ax2.set_title('Anthropogenic propane emissions')
ax2.set_ylabel('Tg yr$^{-1}$')

plt.legend()
plt.savefig( 'plots/emissions.ala-dalsoren.CEDS-2021.png')
plt.close()
