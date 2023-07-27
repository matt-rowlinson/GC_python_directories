


import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use("agg")
plt.style.use('seaborn-darkgrid')

f=plt.figure(figsize=(15,4))
ax1 = f.add_subplot(1,2,1)
# Ethane emissions
xticks=['Global','USA','Europe','China','NH','SH','0-30N','30-60N','60-90N', 'India']
ceds  = [  7.06,  0.33, 0.13+0.37,2.67,  5.67, 1.39, 4.11, 1.89,    0.05,    4.11    ]
me    = [ 10.96,  2.06, 0.31+0.86,4.26,  9.57, 1.39, 5.0,  5.31,    0.10,    5.00    ]

ceds21= [ 4.94,   0.25, 0.13+0.34,2.32, 4.12, 0.82, 2.71, 1.74,     0.05,    2.71    ]
me2   = [ 8.97 ,  1.78, 0.29+0.74,4.33, 8.14, 0.82, 3.82, 5.17,     0.09,    3.82    ]

cs = ['#b2df8a','#33a02c','#cab2d6','#6a3d9a','#1b9e77','#d95f02','#7570b3','#e7298a']
ax1.plot( [xticks[0],xticks[0]], [3.36, 10.],marker='_', solid_capstyle='butt', lw=6,c='y', alpha=.3, label='Bottom up estimates')
ax1.scatter( xticks, ceds  , s=82, c=cs[0], marker='o', label='CEDS 2020-08' )
ax1.scatter( xticks, ceds21, s=82, c=cs[1], marker='o', label='CEDS 2021-06' )

ax1.scatter( xticks, me    , s=82, c=cs[2], marker='*', label='Old with scaling')
ax1.scatter( xticks, me2   , s=82, c=cs[3], marker='*', label='New with scaling' )
ax1.set_ylabel('Tg yr$^{-1}$')
ax1.set_title('Anthropogenic ethane emissions')
plt.legend()

ax2 = f.add_subplot(1,2,2)
## Propane emissions
ceds  = [ 6.40, 0.36, 0.19+0.32, 1.90, 5.03, 1.37, 3.47, 1.88, 0.10, 3.47  ]
me    = [ 8.10, 1.22, 0.37+0.63, 2.30, 6.73, 1.37, 3.84, 3.39, 0.17, 3.84]

ceds21= [ 4.13, 0.33, 0.18+0.28, 1.55, 3.51, 0.62, 2.17, 1.63, 0.10, 2.17 ]
me2   = [ 6.49, 1.71, 0.37+0.57, 2.10, 5.87, 0.62, 2.82, 3.69, 0.18, 2.82 ]

ax2.plot( [xticks[0],xticks[0]], [3., 10.],marker='_', solid_capstyle='butt', lw=6,c='y', alpha=.3, label='Bottom up estimates')
ax2.scatter( xticks, ceds , s=82, c=cs[0], marker='o', label='CEDS 2020-08' )
ax2.scatter( xticks, ceds21,s=82, c=cs[1], marker='o', label='CEDS 2021-06' )

ax2.scatter( xticks, me   , s=82, c=cs[2], marker='*', label='Old with scaling')
ax2.scatter( xticks, me2  , s=82, c=cs[3], marker='*', label='New with scaling')
ax2.set_title('Anthropogenic propane emissions')
ax2.set_ylabel('Tg yr$^{-1}$')

plt.legend(fancybox=True, facecolor='w')
plt.savefig( 'plots/CEDS.emissions.wIndia.png')
plt.close()
