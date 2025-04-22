#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

## Define equation for fit
def LinFit(x, a, b):
    return x*a + b    

f = np.loadtxt('Data_Nup116.dat', dtype=object)
names = f[:,0]
DelGhyd = f[:,1].astype(float)
P116 = f[:,2].astype(float)
P98A = f[:,3].astype(float)
logP116 = np.log(P116)
logP98A = np.log(P98A)

## Fit the data to a log relationship
popt116, pcov116 = curve_fit(LinFit, DelGhyd, logP116)
popt98A, pcov98A = curve_fit(LinFit, DelGhyd, logP98A)

## Calculate variance and R^2
pred116 = LinFit(DelGhyd, *popt116)
pred98A = LinFit(DelGhyd, *popt98A)
residual_116 = pred116 - logP116
residual_98A = pred98A - logP98A
ss_res_116 = np.sum(residual_116**2)
ss_res_98A = np.sum(residual_98A**2)
ss_tot_116 = np.sum((logP116 - np.mean(logP116))**2)
ss_tot_98A = np.sum((logP98A - np.mean(logP98A))**2)
R2_116 = 1 - (ss_res_116 / ss_tot_116)
R2_98A = 1 - (ss_res_98A / ss_tot_98A)

## Set up the colorlist to show the residual magnitude
max_res = np.max([np.max(residual_116),np.max(residual_98A)])
min_res = np.min([np.min(residual_116),np.min(residual_98A)])
plus_minus_dist = 2.5
cmap_idx_116 = (residual_116 + plus_minus_dist) / (2*plus_minus_dist)
cmap_idx_98A = (residual_98A + plus_minus_dist) / (2*plus_minus_dist)
clist_116 = plt.cm.seismic(cmap_idx_116)
clist_98A = plt.cm.seismic(cmap_idx_98A)

## Calculate Gtrans from partition coefficients
DelGtrans_116 = -2.5 * logP116
DelGtrans_98A = -2.5 * logP98A

## Calculate solvation free energy 
DelGsolv_116 = DelGhyd + DelGtrans_116
DelGsolv_98A = DelGhyd + DelGtrans_98A

## Set up the plot
fig, ax = plt.subplots(nrows=2, dpi=150, figsize=[3,5])
ax[0].scatter(DelGhyd, DelGtrans_116, s=10, c='g', lw=0.5, edgecolors='k')
ax[0].set_xlabel('')
ax[0].set_xlim(-6550,-4050)
ax[0].set_xticks([-6550,-5715,-4883,-4050])
ax[0].set_xticklabels(['','','',''])
ax[0].set_ylabel(r'$\Delta G_{trans}$ (kJ/mol)', fontsize=10)
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize=9)
#ax[0].set_ylim(0.05, 500)
ax[0].tick_params('both', direction='in')
ax[0].legend(fontsize=7)
ax[1].scatter(DelGhyd, DelGtrans_98A, s=10, c='darkorchid', lw=0.5, edgecolor='k')
ax[1].set_xlabel(r'$\Delta G_{hyd}$ (kJ/mol)', fontsize=10)
ax[1].set_xlim(-6550,-4050)
ax[1].set_xticklabels(ax[1].get_xticks().astype(int), fontsize=9)
ax[1].set_ylabel(r'$\Delta G_{trans}$ (kJ/mol)', fontsize=10)
ax[1].set_yticklabels(ax[1].get_yticks(), fontsize=9)
#ax[1].set_ylim(0.05, 500)
ax[1].tick_params('both', direction='in')
ax[1].legend(fontsize=7)

plt.tight_layout()

## Save the figure
plt.savefig('Gtrans_Ghyd_plots.pdf')
plt.savefig('Gtrans_Ghyd_plots.png')
plt.show()

## Plot bars
fig, ax = plt.subplots(nrows=2, dpi=150, figsize=[5,5])
ax[0].bar(np.arange(17)-0.3, DelGhyd, color='cyan', width=0.3, label=r'$\Delta G_{hyd}$')
ax[0].bar(np.arange(17), DelGsolv_116, color='darkblue', width=0.3, label=r'$\Delta G_{solv}$')
ax[0].bar(np.arange(17)+0.3, DelGtrans_116, color='limegreen', width=0.3, label=r'$\Delta G_{trans}$')
#ax[0].plot([-0.9,16.9,16.9,-0.9,-0.9], [-15,-15,15,15,-15], c='k', ls='--', lw=1)
ax[0].set_xlim(-1, 17)
ax[0].set_xlabel('')
ax[0].set_xticks(np.arange(17))
ax[0].set_xticklabels(17*[''])
#ax[0].set_ylim(-8300, 1800)
ax[0].set_ylim(-15,15)
ax[0].set_ylabel(r'$\Delta$G (kJ/mol)', fontsize=10)
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize=9)
ax[0].tick_params('both', direction='in')
ax[1].bar(np.arange(17)-0.3, DelGhyd, color='cyan', width=0.3)
ax[1].bar(np.arange(17), DelGsolv_98A, color='darkblue', width=0.3)
ax[1].bar(np.arange(17)+0.3, DelGtrans_98A, color='limegreen', width=0.3)
#ax[1].plot([0.9,16.9,16.9,-0.9,-0.9], [-15,-15,15,15,-15], c='k', ls='--', lw=1)
ax[1].set_xlim(-1, 17)
ax[1].set_xlabel('GFP Variant', fontsize=10)
ax[1].set_xticks(np.arange(17))
ax[1].set_xticklabels(names, fontsize=9)
#ax[1].set_ylim(-8300, 1800)
ax[1].set_ylim(-15,15)
ax[1].set_ylabel(r'$\Delta$G (kJ/mol)', fontsize=10)
ax[1].set_yticklabels(ax[0].get_yticks(), fontsize=9)
ax[1].tick_params('both', direction='in')
ax[0].legend(fontsize=8, ncol=3)

plt.rcParams['pdf.fonttype'] = 42
plt.tight_layout()
plt.savefig('Bar_chart.pdf')
plt.savefig('Bar_chart.png')
plt.show()




