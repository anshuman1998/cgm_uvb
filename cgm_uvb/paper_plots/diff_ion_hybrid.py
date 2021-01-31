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

out_fig_name = 'histogram_hybrid.pdf'
figure_size = [14, 4]
fig, (ax1, ax2, ax3)  = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.15)

path  = '/home/vikram/cgm_uvb/cgm_uvb/paper_plots/more_models'

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

den = [1e-3, 1e-4, 1e-5]
num_list = [4, 5, 6]

for num, ax in zip(num_list, (ax1, ax2, ax3)):
    lname = r'True log n$_{\rm H}$'
    color_list = ['dodgerblue', 'orange', 'green']

    for nH, col in zip(den, color_list):
        d = tab.Table.read(path + '/diff_res_{}ions_newhyb.txt'.format(num), format='ascii')

        n2_all_new = []
        z2_all_new = []
        for m, n, in zip(d['Q'][d['true_nH'] == nH], d['Max_diff_nH'][d['true_nH'] == nH]):
            for uvb in all_uvb:
                if uvb == m:
                    n2_all_new.append(n)
                    n2_all.append(n)
        print(len(n2_all))

        median_lab = r'$\Delta_{\rm max}$ (log n$_{H}) =$' + '{:.2f}'.format(np.median(n2_all_new))
        if ax == ax1:
            binwidth = 0.05
            ax.set_ylabel('Frequency')

            p = ax.hist(n2_all_new, bins=np.arange(min(n2_all_new), max(n2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.85,
                        histtype='step', linewidth=3, density=1,
                        label=lname + ' = {:.1f}'.format(np.log10(nH)) ,
                        color=col)

        else:
            binwidth = 0.022
            #ax.set_ylim(0.01, 25)

            p = ax.hist(n2_all_new, bins=np.arange(min(n2_all_new), max(n2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.85,
                        histtype='step', linewidth=3, density=1,
                        color=col)

        #ax.set_ylim(0.001, 15)

        ax.axvline(np.median(n2_all_new), linestyle='--', color=col, linewidth=2, label = median_lab)
        print(np.median(n2_all_new))


        if ax == ax1:
            ax.annotate(r'Using {} ions'.format(num), xy=(0.55, 0.5), xycoords='axes fraction', fontsize=12)
            ax.annotate(r'Hybrid model', xy=(0.55, 0.4), xycoords='axes fraction', fontsize=12)
            ax.annotate(r'(T = 10$^{5.5}$ K)', xy=(0.55, 0.35), xycoords='axes fraction', fontsize=12)


        else:
            ax.annotate(r'Using {} ions'.format(num), xy=(0.05, 0.5), xycoords='axes fraction', fontsize=12)

for ax in (ax1, ax2, ax3):
    # ax.axvline(0.7, linestyle = '--', color = 'cyan')
    # ax.axvline (0.52, linestyle = '--', color = 'cyan')

    # ax.annotate (r'Median' , xy=(0.5, 0.9), xycoords='axes fraction', fontsize=11)
    # ax.annotate (r'$\Delta_{\rm max}$ log n$_{\rm H}$ = 0.7' , xy=(0.5, 0.83), xycoords='axes fraction', fontsize=11)

    # ax.annotate (r'$\Delta_{\rm max}$ log Z = 0.52' , xy=(0.5, 0.83), xycoords='axes fraction', fontsize=11)

    print(np.median(n2_all), np.median(z2_all))

    ax.legend(loc='best', fontsize=12)

    # ax1.set_ylabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
    # ax1.set_xlabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
    # ax1.set_xlim (0.5, 1.02)
    # ax1.set_ylim (0.3, 1.02)
    # ax2.set_xlim (0.5, 1.0)
    # ax3.set_xlim (0.3, 0.9)

    ax.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
    # ax3.set_xlabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
    # ax.set_xscale('log')
    # ax.set_yscale('log')


    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    ax.tick_params(direction='in', length=7, width=1.7)
    ax.tick_params(direction='in', which='minor', length=4, width=1.7)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    # decorating the plot
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.7)
        ax.spines[axis].set_color('k')

fig.savefig(out_fig_name, bbox_inches='tight')
