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

out_fig_name = 'scatter_den_hybrid_all_4_ions.pdf'
figure_size = [14, 4]
fig, (ax1, ax2, ax3)  = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.15)



path  = '/home/vikram/cgm_uvb/cgm_uvb/paper_plots/more_models/new_h'

#d  = tab.Table.read(path + '/diff_res_2ions_new.txt', format = 'ascii')
ks_array = ['14', '15', '16', '17', '18', '19', '20']
all_uvb = ks_array + ['FG20', 'P19']
n_H_array = [1e-3, 1e-4, 1e-5]
n2_ks = []
z2_ks = []
n2_all = []
z2_all = []

"""
for num in ion_num:
    d = tab.Table.read(path + '/diff_res_{}ions_new.txt'.format(num), format='ascii')
    print(num,  len(d))
    for m, n, Z, in zip(d['Q'], d['Max_diff_nH'], d['Max_diff_Z']):
        for uvb in all_uvb:
            if uvb == m:
                n2_all.append(n)
                z2_all.append(Z)

print(len(n2_all))
ax.scatter(n2_all, z2_all, alpha = 0.5, color = 'b')
"""

num = 4
d = tab.Table.read(path + '/diff_res_{}ions_newhyb.txt'.format(num), format='ascii')

color_list = ['dodgerblue', 'orange', 'green']

for den, col in zip(n_H_array, color_list):
    n2_all_new = []
    z2_all_new = []
    for m, n, Z, in zip(d['Q'][d['true_nH'] == den], d['Max_diff_nH'][d['true_nH'] == den], d['Max_diff_Z'][d['true_nH']==den]):
        for uvb in all_uvb:
            if uvb == m:
                n2_all_new.append(n)
                z2_all_new.append(Z)
                n2_all.append(n)
                z2_all.append(Z)
    print(len(n2_all))
    label ='true (log Z, log' + r' n$_{\rm H}$' + ') = (-1, {:.1f})'.format(np.log10(den))
    ax1.scatter(n2_all_new, z2_all_new, alpha = 0.4, label =label, s = 16, color = col)
    ax1.scatter( np.median(n2_all_new), np.median(z2_all_new), marker = '*', s = 300, color = 'red', edgecolor = 'k', zorder =10)
    ax1.plot ([np.median(n2_all_new), np.median(n2_all_new), 0], [-0.04, np.median(z2_all_new),  np.median(z2_all_new)],
              linestyle = '--', linewidth  = 1.5, color = col, alpha = 0.7)
    ax1.set_ylim(-0.04, 0.59)


    binwidth = 0.04
    ax2.set_ylabel('Frequency')

    median_lab = r'$\Delta_{\rm max}$ (log n$_{H})$' + '\n = {:.2f}'.format(np.median(n2_all_new))

    p = ax2.hist(n2_all_new, bins=np.arange(min(n2_all_new), max(n2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.85,
                histtype='step', linewidth=3, density=1,
                label='' ,
                color=col)


    p = ax3.hist(z2_all_new, bins=np.arange(min(z2_all_new), max(z2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.85,
                histtype='step', linewidth=3, density=1,
                label='' ,
                color=col)

    ax2.axvline(np.median(n2_all_new), linestyle='--', color=col, linewidth=2, label=median_lab)
    median_lab = r'$\Delta_{\rm max}$ (log Z)' + '\n = {:.2f}'.format(np.median(z2_all_new))
    ax3.axvline(np.median(z2_all_new), linestyle='--', color=col, linewidth=2, label=median_lab)
    print(np.median(n2_all_new), np.median(z2_all_new))

    ax2.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
    ax3.set_xlabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')

    #ax2.set_xlim(0.35, 1.25)

    """
    if den == 1e-3:
        col = 'blue'
    elif den ==1e-4:
        col = 'orange'
    else:
        col = 'green'
    """

    #ax.plot ( linestyle = '--', linewidth  = 1.5)


print(np.median(n2_all), np.median(z2_all))

for ax in(ax1, ax2, ax3):

    if ax ==ax1:
        ax.annotate('Hybrid\n'+ 'absorber' +r' (10$^{5.5}$ K)'+ '\n(using {} ions)'.format(num), xy=(0.05, 0.67),
                    xycoords='axes fraction', fontsize=12)
        ax.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
        ax.set_ylabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
        ax.set_xlim(0.35, 1.15)
        #ax.set_ylim(0.15, 0.95)
        ax.legend(loc='best', fontsize=9)
    else:
        ax.legend(loc='best', fontsize=9)

        # ax.set_xscale('log')
        # ax.set_yscale('log')

    #if ax == ax2:
    #    ax.set_xlim(0.3, 1)

    # deco
    ax.tick_params(direction='in', length=7, width=1.7)
    ax.tick_params(direction='in', which='minor', length=4, width=1.7)
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    # decorating the plot
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.7)
        ax.spines[axis].set_color('k')

fig.savefig(out_fig_name, bbox_inches='tight')


#import corner
#import numpy as np

#labels = [r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)', r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)']

#samples = (np.vstack([n2_all_new, z2_all_new])).T

#figure = corner.corner(samples,  labels=labels, quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12})
#figure.savefig('test_s.pdf', bbox_inches='tight')


