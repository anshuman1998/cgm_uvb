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
import os


def plot_plots(den=1e-4, met= -1, ion_name  = 'O+6', outpath = '/home/vikram/cloudy_run/more_fig/'):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5

    final_path  = outpath + '/{}'.format(ion_name)
    if not os.path.isdir(final_path):
        os.mkdir(final_path)

    out_fig_name = final_path + '/w_scatter_hybrid_bimodal_with_Ne8_newfile_logZ{}_logN{}.jpg'.format(met, np.log10(den))
    figure_size = [14, 4]
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.15)

    path = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'
    d = tab.Table.read(path + '/all_combined_logT550.fits')


    ion_num = [5, 8]
    color_list1 = ['dodgerblue', 'green']
    color_list2 = ['orange', 'magenta']
    # , 'grey', 'green', 'magenta',  'cyan', 'gold']

    left_hand = False
    num = ion_num[0]
    ions_used_coloum = d['ions']
    flag = []
    for i in range(len(d)):
        ions = ions_used_coloum[i].split('_')
        ne8 = False
        for j in ions:
            #if j == 'Ne+7':
            if j == ion_name:
                ne8 = True
        if ne8 == True:
            flag.append(1)
        else:
            flag.append(0)

    flag = np.array(flag)
    for num, colr1, colr2, in zip(ion_num, color_list1, color_list2):
        dlow = d[flag == 1]
        dhigh = d[flag == 0]

        # low
        dnum = dlow[dlow['n_ions'] == num]
        sort_d = dnum[dnum['true_nH'] == den]
        sort_d = sort_d[sort_d['true_logZ'] == met]

        n2_all_new = sort_d['n_dmax_KS']
        z2_all_new = sort_d['z_dmax_KS']

        # high
        dnum = dhigh[dhigh['n_ions'] == num]
        sort_d = dnum[dnum['true_nH'] == den]
        sort_d = sort_d[sort_d['true_logZ'] == met]

        n2_all_new_high = sort_d['n_dmax_KS']
        z2_all_new_high = sort_d['z_dmax_KS']

        # ------
        binwidth = 0.02
        binwidth_z = 0.004

        ax1.scatter(n2_all_new, z2_all_new, alpha=0.5, label='{} ions with {}'.format(num, ion_name), s=13, color=colr1)

        ax2.hist(n2_all_new, bins=np.arange(np.min(n2_all_new), np.max(n2_all_new) + binwidth, binwidth), alpha=0.75,
                 histtype='bar', edgecolor='k', linewidth=1.5, density=1, label='{} ions with {}'.format(num, ion_name), color=colr1)
        ax3.hist(z2_all_new, bins=np.arange(np.min(z2_all_new), np.max(z2_all_new) + binwidth_z, binwidth_z), alpha=0.8,
                 histtype='bar', edgecolor='k', linewidth=1.5, density=1, label='{} ions with {}'.format(num, ion_name), color=colr1)

        ax1.scatter(n2_all_new_high, z2_all_new_high, alpha=0.5, label='', s=13, color=colr2)

        ax2.hist(n2_all_new_high, bins=np.arange(np.min(n2_all_new_high), np.max(n2_all_new_high) + binwidth, binwidth),
                 alpha=0.75,
                 histtype='bar', edgecolor='k', linewidth=1.5, density=1, label='', color=colr2)
        ax3.hist(z2_all_new_high, bins=np.arange(np.min(z2_all_new_high), np.max(z2_all_new_high) + binwidth_z, binwidth_z),
                 alpha=0.8,
                 histtype='bar', edgecolor='k', linewidth=1.5, density=1, label='', color=colr2)

        print(np.median(n2_all_new), np.median(z2_all_new), 'for ions', num)

    ax1.annotate('Hybrid \n' + r'absorber (10$^{5.5}$ K) ', xy=(0.06, 0.75), xycoords='axes fraction', fontsize=12)
    ax1.annotate(r'True (log Z, log n$_{\rm H}$) = ' + '({:.0f}, {:.0f})'.format(met, np.log10(den)), xy=(0.06, 0.67),
                 xycoords='axes fraction', fontsize=12)

    ax2.axvline(np.median(n2_all_new), linestyle='--', color='blue', zorder=20)
    ax3.axvline(np.median(z2_all_new), linestyle='--', color='blue', zorder=20)
    ax2.axvline(np.median(n2_all_new_high), linestyle='--', color='red', zorder=20)
    ax3.axvline(np.median(z2_all_new_high), linestyle='--', color='red', zorder=20)

    ax2.annotate(r'Median', xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
    ax2.annotate(r'$\Delta_{\rm max}$ (log n$_{\rm H}$)' + '= {:.2f}'.format(np.median(n2_all_new)), xy=(0.6, 0.83),
                 xycoords='axes fraction', fontsize=11)
    ax2.annotate(r'$\Delta_{\rm max}$ (log n$_{\rm H}$)' + '= {:.2f}'.format(np.median(n2_all_new_high)),
                 xy=(0.6, 0.73), xycoords='axes fraction', fontsize=11)

    ax3.annotate(r'Median', xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
    ax3.annotate(r'$\Delta_{\rm max}$' + ' (log Z) = {:.2f}'.format(np.median(z2_all_new)), xy=(0.6, 0.83),
                 xycoords='axes fraction', fontsize=11)
    ax3.annotate(r'$\Delta_{\rm max}$' + ' (log Z) = {:.2f}'.format(np.median(z2_all_new_high)), xy=(0.6, 0.73),
                 xycoords='axes fraction', fontsize=11)

    ##print(np.median(n2_all), np.median(z2_all))

    ax1.legend(loc='best', fontsize=12, ncol=2)

    ax1.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
    ax1.set_ylabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
    ax1.set_xlim (-0.1, 1.1)
    ax1.set_ylim (-0.1, 0.6)
    # ax2.set_xlim (0.45, 1.07)
    #ax3.set_xlim(-0.05, 0.55)

    ax2.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
    ax3.set_xlabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
    ax2.set_ylabel('Frequency')
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

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
    plt.close()

#plt.show()

Zarray = [-2, -1, 0]
Narray = [1e-5, 1e-4, 1e-3]
ion_names = [ 'C+2','C+3','N+2','N+4','O+2','O+3','O+4','Si+3', 'O+5', 'Ne+7', 'N+2', 'N+3', 'N+4', 'O+2', 'S+4']
for ion in ion_names:
    for n in Narray:
        for z in Zarray:
            plot_plots(den=n, met=z, ion_name= ion)