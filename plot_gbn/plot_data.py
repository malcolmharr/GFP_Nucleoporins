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
"""pred116 = np.polyfit(DelGhyd, logP116, 1)
pred98A = np.polyfit(DelGhyd, logP98A, 1)
func116 = np.poly1d(pred116)
func98A = np.poly1d(pred98A)"""
popt116, pcov116 = curve_fit(LinFit, DelGhyd, logP116)
popt98A, pcov98A = curve_fit(LinFit, DelGhyd, logP98A)


## Set up the plot
fig, ax = plt.subplots(nrows=2, dpi=300, figsize=[3,5])
ax[0].scatter(DelGhyd, P116, s=10, c='g', lw=0.5, edgecolors='k')
ax[0].plot([-4650,-4050], np.exp(LinFit(np.array([-4650,-4050]), *popt116)), c='g', lw=1, label=r'ln(P)=%.4f$\Delta G_{hyd}$+%.1f'%(popt116[0], popt116[1]))
ax[0].set_xlabel('')
ax[0].set_xlim(-4650,-4050)
ax[0].set_xticks([-4650,-4450,-4250,-4050])
ax[0].set_xticklabels(['','','',''])
ax[0].set_ylabel(r'P$_{water->Nup116}$', fontsize=10)
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize=9)
ax[0].set_yscale('log')
ax[0].set_ylim(0.05, 500)
ax[0].tick_params('both', direction='in')
ax[0].legend(fontsize=7)
ax[1].scatter(DelGhyd, P98A, s=10, c='darkorchid', lw=0.5, edgecolors='k')
ax[1].plot([-4650,-4050], np.exp(LinFit(np.array([-4650,-4050]), *popt98A)), c='darkorchid', lw=1, label=r'ln(P)=%.4f$\Delta G_{hyd}$+%.1f'%(popt98A[0], popt98A[1]))
ax[1].set_xlabel(r'$\Delta G_{hyd}$ (kJ/mol)', fontsize=10)
ax[1].set_xlim(-4650,-4050)
ax[1].set_xticklabels(ax[1].get_xticks().astype(int), fontsize=9)
ax[1].set_ylabel(r'P$_{water->MacNup98A}$', fontsize=10)
ax[1].set_yticklabels(ax[1].get_yticks(), fontsize=9)
ax[1].set_yscale('log')
ax[1].set_ylim(0.05, 500)
ax[1].tick_params('both', direction='in')
ax[1].legend(fontsize=7)

plt.tight_layout()

plt.rcParams['pdf.fonttype'] = 42

## Save the figure
plt.savefig('P_Ghyd_plots_test.pdf')
plt.savefig('P_Ghyd_plots_test.png')
plt.show()
