
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use("agg")
plt.style.use('seaborn-darkgrid')

f=plt.figure(figsize=(10,5))
ax1 = f.add_subplot(1,2,1)
# Ethane emissions
xticks=['Global','USA','Europe','China']
ceds  = [  4.94,  0.25, 0.13+0.37, 2.32 ]
edgar = [ 4.855, 0.29, 0.50, 0.65 ]
me   = [ 8.97 , 1.78, 0.29+0.74, 4.33 ]
dalso = [ 13.60, np.nan, np.nan, np.nan ]

cs = ['#1b9e77','#d95f02','#7570b3','#e7298a']
ax1.plot( [xticks[0],xticks[0]], [3.36, 10.],marker='_', solid_capstyle='butt', lw=6,c='y', alpha=.3, label='Bottom up estimates')
ax1.scatter( xticks, ceds , s=42, c=cs[0], marker='o', label='CEDS' )
ax1.scatter( xticks, edgar, s=42, c=cs[1], marker='^', label='EGDAR')
ax1.scatter( xticks, dalso, s=52, c=cs[3], marker='+', label='Dalsøren et al. (2019)')
ax1.scatter( xticks, me   , s=82, c=cs[2], marker='*', label='This study')
ax1.set_ylabel('Tg yr$^{-1}$')
ax1.set_title('Anthropogenic ethane emissions')
plt.legend()

ax2 = f.add_subplot(1,2,2)
## Propane emissions
ceds  = [ 4.13, 0.33, 0.18+0.28, 1.55 ]
edgar = [ 6.22, 0.2379, 0.411, 1.4 ]
me    = [ 6.49, 1.71, 0.37+0.57, 2.10 ]
dalso = [ 16.10, np.nan, np.nan, np.nan ]

ax2.plot( [xticks[0],xticks[0]], [3., 10.],marker='_', solid_capstyle='butt', lw=6,c='y', alpha=.3, label='Bottom up estimates')
ax2.scatter( xticks, ceds , s=42, c=cs[0], marker='o', label='CEDS' )
ax2.scatter( xticks, edgar, s=42, c=cs[1], marker='^', label='EDGAR')
ax2.scatter( xticks, dalso, s=52, c=cs[3], marker='+', label='Dalsøren et al. (2019)')
ax2.scatter( xticks, me   , s=82, c=cs[2], marker='*', label='This study')
ax2.set_title('Anthropogenic propane emissions')
ax2.set_ylabel('Tg yr$^{-1}$')

plt.tight_layout()
plt.legend()
plt.savefig( 'plots/figure_2.png')
plt.close()
