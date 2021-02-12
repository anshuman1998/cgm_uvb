import astropy.table as t
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
from matplotlib.ticker import ScalarFormatter
from astropy.cosmology import Planck15 as planck
import scipy.optimize as sciopt
import astropy.table as tab


# setting the figure
font = {'family': 'serif', 'weight': 'normal', 'size': 14}
plt.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5

out_fig_name = 'scatter_hybrid_new.pdf'
figure_size = [14, 4]
fig, (ax1, ax2, ax3)  = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.15)


path  = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'


den = 1e-4
met= -1
ion_num = [3, 5,  8]

d = tab.Table.read(path + '/all_combined_logT550.fits')

for num in ion_num :
    dnum = d[d['n_ions'] == num]
    sort_d = dnum[dnum['true_nH'] == den]
    sort_d = sort_d[sort_d['true_logZ'] == met]

    n2_all_new = sort_d['n_dmax']
    z2_all_new = sort_d['z_dmax']

    ax1.scatter(n2_all_new, z2_all_new, alpha=0.5, label='{} ions'.format(num), s=13)

    binwidth = 0.035
    ax2.hist(n2_all_new, bins=np.arange(min(n2_all_new), max(n2_all_new) + binwidth, binwidth), alpha = 0.75,
             histtype = 'bar', edgecolor= 'k', linewidth = 1.5, density = 1, label = '{} ions'.format(num) )
    binwidth = 0.025
    ax3.hist(z2_all_new, bins=np.arange(min(z2_all_new), max(z2_all_new) + binwidth, binwidth), alpha = 0.8,
             histtype = 'bar',  edgecolor= 'k', linewidth = 1.5, density = 1, label = '{} ions'.format(num))

    print(np.median(n2_all_new), np.median(z2_all_new), 'for ions', num)

ax1.annotate ('Hybrid \n' + r'absorber (10$^{5.5}$ K) ' , xy=(0.06, 0.75), xycoords='axes fraction', fontsize=12)

ax1.annotate (r'True (log Z, log n$_{\rm H}$) = (-1, -4)' , xy=(0.06, 0.67), xycoords='axes fraction', fontsize=12)


ax2.axvline(np.median(n2_all_new), linestyle = '--', color = 'cyan')
ax3.axvline (np.median(z2_all_new), linestyle = '--', color = 'cyan')

ax2.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
ax2.annotate (r'$\Delta_{\rm max}$ (log n$_{\rm H}$)' + '= {:.2f}'.format(np.median(n2_all_new)) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)


ax3.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
ax3.annotate (r'$\Delta_{\rm max}$' +' (log Z) = {:.2f}'.format(np.median(z2_all_new)) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)


print(np.median(n2_all_new), np.median(z2_all_new))

ax1.legend(loc = 'best',  fontsize = 12, ncol=2)

ax1.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
ax1.set_ylabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
#ax1.set_xlim (0.5, 1.4)
#ax1.set_ylim (0.25, 1.02)
#ax2.set_xlim (0.45, 1.4)
#ax3.set_xlim (0.25, 1.2)

ax2.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
ax3.set_xlabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
ax2.set_ylabel('Frequency')
#ax.set_xscale('log')
#ax.set_yscale('log')


#ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

for ax in (ax1, ax2, ax3):
    # deco
    ax.tick_params(direction='in', length=7, width=1.7)
    ax.tick_params(direction='in', which='minor', length=4, width=1.7)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    # decorating the plot
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.7)
        ax.spines[axis].set_color('k')

fig.savefig(out_fig_name, bbox_inches='tight')

#plt.show()