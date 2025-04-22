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
print(f"R^2 value for 116: {R2_116}")
print(f"R^2 value for 98A: {R2_98A}")

## Set up the colorlist to show the residual magnitude
max_res = np.max([np.max(residual_116),np.max(residual_98A)])
min_res = np.min([np.min(residual_116),np.min(residual_98A)])
plus_minus_dist = 2.5
cmap_idx_116 = (residual_116 + plus_minus_dist) / (2*plus_minus_dist)
cmap_idx_98A = (residual_98A + plus_minus_dist) / (2*plus_minus_dist)
clist_116 = plt.cm.jet_r(cmap_idx_116)
clist_98A = plt.cm.jet_r(cmap_idx_98A)

## Set up the plot
fig, ax = plt.subplots(nrows=2, dpi=300, figsize=[3,4])
ax[0].scatter(DelGhyd, logP116, s=10, c='g', lw=0.5, edgecolors='k')
ax[0].plot([-8500,-7200], LinFit(np.array([-8500,-7200]), *popt116), c='g', lw=1, label=r'ln(P)=%.4f$\Delta G_{hyd}$+%.1f'%(popt116[0], popt116[1]))
ax[0].plot([-8500,-7200], np.exp((np.array([-8500,-7200])/2.5)), c='k', lw=0.75, ls='--')
## Plot residuals
for i in range(len(residual_116)):
    print(i)
    ax[0].plot([DelGhyd[i]]*2, [pred116[i],logP116[i]], c=clist_116[i], lw=0.5)
ax[0].set_xlabel('')
ax[0].set_xlim(-8500,-7200)
ax[0].set_xticks([-8400,-8000,-7600,-7200])
ax[0].set_xticklabels(['','','',''])
ax[0].set_ylabel(r'ln(P$_{water->Nup116}$)', fontsize=10)
ax[0].set_yticklabels(ax[0].get_yticks(), fontsize=9)
ax[0].tick_params('both', direction='in')
#ax[0].legend(fontsize=7)
ax[1].scatter(DelGhyd, logP98A, s=10, c='darkorchid', lw=0.5, edgecolor='k')
ax[1].plot([-8500,-7200], LinFit(np.array([-8500,-7200]), *popt98A), c='darkorchid', lw=1, label=r'ln(P)=%.4f$\Delta G_{hyd}$+%.1f'%(popt98A[0], popt98A[1]))
## Plot residuals
for i in range(len(residual_116)):
    print(i)
    ax[1].plot([DelGhyd[i]]*2, [pred98A[i],logP98A[i]], c=clist_98A[i], lw=0.5)
ax[1].set_xlabel(r'$\Delta G_{hyd}$ (kJ/mol)', fontsize=10)
ax[1].set_xlim(-8500,-7200)
ax[1].set_xticklabels(ax[1].get_xticks().astype(int), fontsize=9)
ax[1].set_ylabel(r'ln(P$_{water->MacNup98A}$)', fontsize=10)
ax[1].set_yticklabels(ax[1].get_yticks(), fontsize=9)
ax[1].tick_params('both', direction='in')
#ax[1].legend(fontsize=7)

plt.tight_layout()
plt.rcParams['pdf.fonttype'] = 42
## Save the figure
plt.savefig('P_Ghyd_plots.pdf')
plt.savefig('P_Ghyd_plots.png')
plt.show()

## Save residuals into new file
outp = np.zeros((len(residual_116), 3), dtype=object)
outp[:,0] = names
outp[:,1] = residual_116
outp[:,2] = residual_98A
np.savetxt('Residual_values.dat', outp, fmt='%8s %8.4f %8.4f', header='Name    116_res  98A_res')
R2_out = np.array([R2_116, R2_98A]).reshape(1, -1)  # Reshaping to a 1x2 array
np.savetxt('Coefficient_of_determination.txt', R2_out, fmt='%8.4f %8.4f', header='R2_116  R2_98A')
