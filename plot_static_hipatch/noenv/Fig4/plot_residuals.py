#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

## Import the data
f = np.loadtxt('Residual_values.dat', dtype='object')
names = f[:,0]
res_116 = f[:,1].astype(float)
res_98A = f[:,2].astype(float)

## Sort the data by residual value
names_116 = names[res_116.argsort()]
res_116 = res_116[res_116.argsort()]
names_98A = names[res_98A.argsort()]
res_98A = res_98A[res_98A.argsort()]

## Define colorlists same as previous script
plus_minus_dist = 2.5
cmap_idx_116 = (res_116 + plus_minus_dist) / (2*plus_minus_dist)
cmap_idx_98A = (res_98A + plus_minus_dist) / (2*plus_minus_dist)
clist_116 = plt.cm.jet_r(cmap_idx_116)
clist_98A = plt.cm.jet_r(cmap_idx_98A)

## Plot the data
fig, ax = plt.subplots(nrows=2, dpi=150, figsize=[5,5])
ax[0].bar(np.arange(17), res_116, width=0.9, color=clist_116)
ax[0].set_xlim(-0.6,16.6)
ax[0].set_xlabel('')
ax[0].set_xticks(np.arange(17))
ax[0].set_xticklabels(names_116, fontsize=9)
ax[0].set_ylabel(r'$\Delta$ ln(P)', fontsize=10)
ax[0].set_yticks([-3,0,3])
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize=9)
ax[0].tick_params('both', direction='in')
ax[1].bar(np.arange(17), res_98A, width=0.9, color=clist_98A)
ax[1].set_xlim(-0.6,16.6)
ax[1].set_xlabel('GFP variant', fontsize=10)
ax[1].set_xticks(np.arange(17))
ax[1].set_xticklabels(names_98A, fontsize=9)
ax[1].set_ylabel(r'$\Delta$ ln(P)', fontsize=10)
ax[1].set_yticks([-3,0,3])
ax[1].set_yticklabels(ax[1].get_yticks(), fontsize=9)
ax[1].tick_params('both', direction='in')

plt.tight_layout()
plt.savefig('Residuals.pdf')
plt.savefig('Residuals.png')
plt.show()

