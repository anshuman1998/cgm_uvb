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

def make_many_plots(den, met, outpath, ion_num = [5], uvb_true ='KS18', q_true = 18):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5


    final_path  = outpath
    if not os.path.isdir(final_path):
        os.mkdir(final_path)

    out_fig_name = outpath + '/pg116.pdf'.format(uvb_true, q_true, ion_num[0], met, np.log10(den))
    #figure_size = [14, 4]
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))

    figure_size = [5, 10]
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, figsize=(figure_size[0], figure_size[1]))


    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0.2, wspace=0.0)

    path = '/home/vikram/cloudy_run/diff_op/photo_NH15'

    #uvb = ['KS18', 'HM12', 'P19', 'FG20']
    uvb = ['KS18']

    uvb_Q = [14, 15, 16, 17, 18, 19, 20]

    #uvb_Q = [14, 20]

    uvb_models = []
    the_Q_values = []
    for background in uvb:
        if background == 'KS18':
            for q in uvb_Q:
                uvb_models.append(background)
                the_Q_values.append(q)
        else:
            q = 18
            uvb_models.append(background)
            the_Q_values.append(q)

    # ks_array = ['14', '15', '16', '17', '18', '19', '20']
    # all_uvb = ks_array + ['FG20', 'P19']
    # n_H_array = [1e-3, 1e-4, 1e-5]

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

    binwidth = 0.2



    color_list = ['dodgerblue', 'orange', 'green', 'magenta']

    num=ion_num[0]

    for the_uvb, q_val in zip(uvb_models, the_Q_values):



        output_file = '/home/vikram/cgm_uvb/hst/pg116_Q{}.fits'.format(q_val)
        sort_d = tab.Table.read(output_file)
        n_max_array = sort_d['nH']
        z_max_array = sort_d['Z']






        if the_uvb =='KS18':
            label0 = 'KS19 Q{}'.format(the_uvb, q_val)
        else:
            label0 = '{} '.format(the_uvb)

        sc= ax1.scatter(n_max_array, z_max_array, alpha=0.35, s=13)

        plt.draw()
        color_draw = sc.get_edgecolor()[0]
        new_colr = [color_draw[0], color_draw[1], color_draw[2], 1]
        #print(line[0].get_color())
        new_colr_alp = [color_draw[0], color_draw[1], color_draw[2], 0.7]
        ax1.scatter(n_max_array, z_max_array, alpha=0.3, facecolor = 'None', s=13, edgecolor  = new_colr_alp, linewidth = 0.6, zorder = -1)

        ax2.hist(n_max_array, bins=np.arange(min(n_max_array)-0.5*binwidth, max(n_max_array) + 1.5*binwidth, binwidth),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = 'k')

        ax3.hist(z_max_array, bins=np.arange(min(z_max_array)-0.5*binwidth, max(z_max_array) + 1.5*binwidth, binwidth),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = 'k')

        n_med = np.median(n_max_array)
        z_med = np.median(z_max_array)

        label_name = r'log n$_{\rm H}$) =' + '{:.2f} ->{} Q{}'.format(n_med, the_uvb, q_val)
        #ax2.axvline(n_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(n_med)
        if the_uvb != 'P19':
            ax2.annotate(txt, (n_med+0.008, 90), fontsize=8, color = new_colr)
        else:
            ax2.annotate(txt, (n_med+0.008, 90), fontsize=8, color = new_colr)

        label_name = r' log Z =' + '{:.2f} ->{} Q{}'.format(z_med, the_uvb, q_val)
        #ax3.axvline(z_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(z_med)
        if the_uvb != 'P19' and q_val !=20:
            ax3.annotate(txt, (z_med+0.008, 150), fontsize=8, color = new_colr)
        else:
            ax3.annotate(txt, (z_med+0.008, 135), fontsize=8, color = new_colr)







        print(the_uvb)



    #ax1.annotate('Photoionized absorber \n' +'True UVB KS19 (Q18)', xy=(0.06, 0.83), xycoords='axes fraction', fontsize=11)

    #ax1.annotate(r'True (log Z, log n$_{\rm H}$) = ' + '({:.0f}, {:.0f})'.format(met, np.log10(den)),
    #             xy=(0.06, 0.76), xycoords='axes fraction', fontsize=11)

    #x = [-1000, -1000]
    #ax1.scatter(x, x, alpha=0.5, s=18, label='Using combinations of 5 ions out of 8', color = 'b')
    #lg = ax1.legend(fontsize=10, loc=4)
    #lg.get_frame().set_facecolor('none')

    """
    n_med = np.median(n_max_array)
    z_med = np.median(z_max_array)

    ax2.axvline(n_med, linestyle = '--', color = 'cyan')
    ax3.axvline (z_med, linestyle = '--', color = 'cyan')

    ax2.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
    ax2.annotate (r'$\Delta_{\rm max}$ (log n$_{\rm H}$) =' + '{:.2f}'.format(n_med) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)


    ax3.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
    ax3.annotate (r'$\Delta_{\rm max}$ (log Z) =' + '{:.2f}'.format(z_med) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)


    print(n_med, z_med)
    """

    #ax1.legend(loc='best', fontsize=9, ncol=2)
    #ax2.legend(loc='best', fontsize=9)
    #ax3.legend(loc='best', fontsize=9)

    #ax1.set_xlabel(r'log n$_{\rm H}$ (cm $^{-3}$)')
    ax1.set_ylabel(r'log Z(Z$_{\odot}$)')

    add_num = 0.49
    #ax1.set_xlim(-4.4, np.log10(den)+add_num)
    #ax1.set_ylim(-1.4, met+add_num)
    #ax2.set_xlim(-4.4, np.log10(den)+add_num)
    #ax3.set_xlim(-1.4, met+add_num)
    #ax2.set_ylim(0, 65)
    #ax3.set_ylim(0, 65)

    ax2.set_xlabel(r'log n$_{\rm H}$ (cm $^{-3}$)')
    ax3.set_xlabel(r'log Z(Z$_{\odot}$)')
    ax2.set_ylabel('Frequency')
    ax3.set_ylabel('Frequency')

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


    return


"""
Zarray = [-2, -1, 0]
Narray = [1e-5, 1e-4, 1e-3]

for n in Narray:
    for z in Zarray:
        make_many_plots(den=n, met=z, outpath= '/home/vikram/cloudy_run/more_fig/photo')
"""

Zarray = [-1]
Narray = [1e-4]

for n in Narray:
    for z in Zarray:
        make_many_plots(den=n, met=z, outpath= '/home/vikram/Work/my_papers/proposals_for_grants_and_time/hst31')