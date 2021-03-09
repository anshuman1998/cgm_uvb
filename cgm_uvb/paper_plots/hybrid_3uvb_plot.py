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

def make_neon_plots(den, met, outpath, ion_num = [8], uvb_true ='KS18', q_true = 18):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5


    final_path  = outpath
    if not os.path.isdir(final_path):
        os.mkdir(final_path)

    out_fig_name = outpath + '/paper_separate_true_UVB_{}_Q{}_{:.0f}ions_logZ{}_logN{}_hybrid_T550.pdf'.format(uvb_true, q_true, ion_num[0], met, np.log10(den))
    #figure_size = [14, 4]
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))

    figure_size = [5, 10]
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, figsize=(figure_size[0], figure_size[1]))


    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0.2, wspace=0.0)

    path = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'

    #uvb = ['KS18', 'HM12', 'P19', 'FG20']
    uvb = ['KS18',  'FG20']

    uvb_Q = [14, 16, 18]

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

    binwidth = 0.02
    binwidth_met = 0.004

    d = tab.Table.read(path + '/all_combined_logT550.fits')
    d = d [d['true_uvb'] == uvb_true]
    d = d [d['true_Q'] == q_true]

    color_list = ['dodgerblue', 'orange', 'green', 'magenta']

    num=ion_num[0]
    dnum = d[d['n_ions'] == num]
    sort_d = dnum[dnum['true_nH'] == den]
    sort_d = sort_d[sort_d['true_logZ'] == met]

    ions_used_coloum = sort_d['ions']
    flag = []
    for i in range(len(sort_d)):
        ions = ions_used_coloum[i].split('_')
        ne8 = False
        for j in ions:
            if j == 'Ne+7':
                ne8 = True
        if ne8 == True:
            flag.append(1)
        else:
            flag.append(0)

    flag = np.array(flag)

    sort_d_with_ne = sort_d[flag==1]
    sort_d_no_ne = sort_d[flag==0]


    for the_uvb, q_val in zip(uvb_models, the_Q_values):

        col_name_nH = 'nH_' + the_uvb + '_Q{}'.format(q_val)
        col_name_Z = 'Z_' + the_uvb + '_Q{}'.format(q_val)


        n_max_array = sort_d_with_ne[col_name_nH]
        z_max_array = sort_d_with_ne[col_name_Z]

        n_max_array_none= sort_d_no_ne[col_name_nH]
        z_max_array_none = sort_d_no_ne[col_name_Z]


        print(len(n_max_array), ': num', len(n_max_array_none))

        if the_uvb =='KS18':
            label0 = 'KS19 Q{}'.format(the_uvb, q_val)
        else:
            label0 = '{} '.format(the_uvb)

        sc= ax1.scatter(n_max_array, z_max_array, alpha=0.3, s=13)

        plt.draw()
        color_draw = sc.get_edgecolor()[0]
        new_colr = [color_draw[0], color_draw[1], color_draw[2], 1]
        #print(line[0].get_color())
        new_colr_alp = [color_draw[0], color_draw[1], color_draw[2], 0.7]
        #--------------------second hist
        ax1.scatter(n_max_array_none, z_max_array_none, s=13, marker= 'D', color= 'cyan', alpha = 0.3)
        if q_val == 18:
            sc = ax1.scatter(n_max_array, z_max_array, alpha=0.3, s=13, facecolor  = new_colr_alp)

        ax1.scatter(n_max_array_none, z_max_array_none, alpha=0.3, facecolor = 'None', s=13, edgecolor  = new_colr_alp,
                    linewidth = 0.3, zorder = -1, marker = 'D')
        ax2.hist(n_max_array_none, bins=np.arange(min(n_max_array_none)-0.5*binwidth, max(n_max_array_none) + 1.5*binwidth, binwidth),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = new_colr_alp, facecolor= 'cyan', hatch = '/')

        ax3.hist(z_max_array_none, bins=np.arange(min(z_max_array_none)-0.5*binwidth, max(z_max_array_none) + 1.5*binwidth_met, binwidth_met),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = new_colr_alp, facecolor = 'cyan', hatch = '/')



        ax1.scatter(n_max_array, z_max_array, alpha=0.3, facecolor = 'None', s=13, edgecolor  = new_colr_alp, linewidth = 0.6, zorder = -1)

        ax2.hist(n_max_array, bins=np.arange(min(n_max_array)-0.5*binwidth, max(n_max_array) + 1.5*binwidth, binwidth),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = new_colr_alp, facecolor = new_colr_alp)

        ax3.hist(z_max_array, bins=np.arange(min(z_max_array)-0.5*binwidth, max(z_max_array) + 1.5*binwidth_met, binwidth_met),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = new_colr_alp, facecolor =new_colr_alp)


        n_med = np.median(n_max_array)
        z_med = np.median(z_max_array)

        n_med_none = np.median(n_max_array_none)
        z_med_none = np.median(z_max_array_none)

        label_name = r'log n$_{\rm H}$) =' + '{:.2f} ->{} Q{}'.format(n_med, the_uvb, q_val)
        ax2.axvline(n_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(n_med)
        if the_uvb != 'P19':
            ax2.annotate(txt, (n_med+0.015, 80), fontsize=8, color = new_colr)
        else:
            ax2.annotate(txt, (n_med+0.015, 135), fontsize=8, color = new_colr)


        txt = '{:.2f}'.format(n_med_none)
        ax2.axvline(n_med_none, label=label_name, linewidth = 0.8, dashes=[5, 7], color = 'cyan', zorder =50)
        if q_val != 14:
            ax2.annotate(txt, (n_med_none+0.018, 50), fontsize=8, color = 'cyan')
        else:
            ax2.annotate(txt, (n_med_none+0.025, 50), fontsize=8, color = 'cyan')

        label_name = r' log Z =' + '{:.2f} ->{} Q{}'.format(z_med, the_uvb, q_val)
        ax3.axvline(z_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(z_med)
        if q_val == 16:
            ax3.annotate(txt, (z_med+0.002, 75), fontsize=8, color = new_colr)
        else:
            ax3.annotate(txt, (z_med+0.002, 80), fontsize=8, color = new_colr)

        ax3.axvline(z_med_none, label=label_name, linewidth = 0.8, dashes=[5, 7], color = 'cyan', zorder =50)
        txt = '{:.2f}'.format(z_med_none)
        if q_val == 16:
            ax3.annotate(txt, (z_med_none+0.002, 55), fontsize=8, color = 'cyan')
        else:
            ax3.annotate(txt, (z_med_none+0.002, 60), fontsize=8, color = 'cyan')


        if the_uvb != 'P19' and the_uvb !='FG20':
            txt = 'Q{}'.format(q_val)
            if q_val <= 16:
                ax1.annotate(txt, (n_med+ 0.05, z_med + 0.022), fontsize=8, color = new_colr)
            else:
                ax1.annotate(txt, (n_med - 0.02, z_med + 0.022), fontsize=8, color=new_colr)
        else:
            txt = '{}'.format(the_uvb)
            if the_uvb == 'FG20':
                ax1.annotate(txt, (n_med -0.03, z_med + 0.05), fontsize=8, color = new_colr)
            else:
                ax1.annotate(txt, (n_med -0.03, z_med - 0.1), fontsize=8, color = new_colr)

        print(the_uvb)



    ax1.annotate(r'Hybrid absorber 10$^{5.5}$K' +'\nTrue UVB KS19 (Q18)', xy=(0.11, 0.83), xycoords='axes fraction', fontsize=11)

    ax1.annotate(r'True (log Z, log n$_{\rm H}$) = ' + '({:.0f}, {:.0f})'.format(met, np.log10(den)),
                 xy=(0.11, 0.76), xycoords='axes fraction', fontsize=11)

    x = [-1000, -1000]
    ax1.scatter(x, x, alpha=0.5, s=18, label='with Ne VIII', color = 'b', facecolor = 'white')

    x = [-1000, -1000]
    ax1.scatter(x, x, alpha=0.5, s=18, marker = 'D', label='without Ne VIII', color = 'cyan')
    lg = ax1.legend(fontsize=11, loc=4)
    lg.get_frame().set_facecolor('none')

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

    add_num0 = 0.42
    ax1.set_xlim(np.log10(den)-add_num0, np.log10(den)+add_num0 + 0.22)
    add_num  = 0.15
    #ax1.set_ylim(-1.09, -0.88)
    ax1.set_ylim(-1.065, -0.885)

    ax2.set_xlim(np.log10(den)-add_num0, np.log10(den)+add_num0 +0.22)
    ax3.set_xlim(-1.05, -0.93)
    #ax2.set_ylim(0, 174)
    #ax3.set_ylim(0, 174)

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


def make_normal_plots(den, met, outpath, ion_num = [8], uvb_true ='KS18', q_true = 18):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5


    final_path  = outpath
    if not os.path.isdir(final_path):
        os.mkdir(final_path)

    out_fig_name = outpath + '/paper_normal_true_UVB_{}_Q{}_{:.0f}ions_logZ{}_logN{}_hybrid_T550.pdf'.format(uvb_true, q_true, ion_num[0], met, np.log10(den))
    #figure_size = [14, 4]
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))

    figure_size = [5, 10]
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, figsize=(figure_size[0], figure_size[1]))


    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0.2, wspace=0.0)

    path = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'

    #uvb = ['KS18', 'HM12', 'P19', 'FG20']
    uvb = ['KS18',  'FG20']

    uvb_Q = [14, 16, 18]

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

    binwidth = 0.02
    binwidth_met = 0.004

    d = tab.Table.read(path + '/all_combined_logT550.fits')
    d = d [d['true_uvb'] == uvb_true]
    d = d [d['true_Q'] == q_true]

    color_list = ['dodgerblue', 'orange', 'green', 'magenta']

    num=ion_num[0]
    dnum = d[d['n_ions'] == num]
    sort_d = dnum[dnum['true_nH'] == den]
    sort_d = sort_d[sort_d['true_logZ'] == met]


    for the_uvb, q_val in zip(uvb_models, the_Q_values):

        col_name_nH = 'nH_' + the_uvb + '_Q{}'.format(q_val)
        col_name_Z = 'Z_' + the_uvb + '_Q{}'.format(q_val)


        n_max_array = sort_d[col_name_nH]
        z_max_array = sort_d[col_name_Z]


        if the_uvb =='KS18':
            label0 = 'KS19 Q{}'.format(the_uvb, q_val)
        else:
            label0 = '{} '.format(the_uvb)

        sc= ax1.scatter(n_max_array, z_max_array, alpha=0.3, s=13)

        plt.draw()
        color_draw = sc.get_edgecolor()[0]
        new_colr = [color_draw[0], color_draw[1], color_draw[2], 1]
        #print(line[0].get_color())
        new_colr_alp = [color_draw[0], color_draw[1], color_draw[2], 0.7]
        #--------------------second hist

        ax1.scatter(n_max_array, z_max_array, alpha=0.3, facecolor = 'None', s=13, edgecolor  = new_colr_alp, linewidth = 0.6, zorder = -1)

        ax2.hist(n_max_array, bins=np.arange(min(n_max_array)-0.5*binwidth, max(n_max_array) + 1.5*binwidth, binwidth),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = new_colr_alp, facecolor = new_colr_alp)

        ax3.hist(z_max_array, bins=np.arange(min(z_max_array)-0.5*binwidth, max(z_max_array) + 1.5*binwidth_met, binwidth_met),
                 alpha=0.75,
                 histtype='stepfilled', linewidth=1.2, edgecolor = new_colr_alp, facecolor =new_colr_alp)


        n_med = np.median(n_max_array)
        z_med = np.median(z_max_array)


        label_name = r'log n$_{\rm H}$) =' + '{:.2f} ->{} Q{}'.format(n_med, the_uvb, q_val)
        ax2.axvline(n_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(n_med)
        if the_uvb != 'P19':
            ax2.annotate(txt, (n_med+0.014, 150), fontsize=8, color = new_colr)
        else:
            ax2.annotate(txt, (n_med+0.014, 135), fontsize=8, color = new_colr)


        label_name = r' log Z =' + '{:.2f} ->{} Q{}'.format(z_med, the_uvb, q_val)
        ax3.axvline(z_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(z_med)
        if q_val == 16:
            ax3.annotate(txt, (z_med+0.002, 150), fontsize=8, color = new_colr)
        else:
            ax3.annotate(txt, (z_med+0.002, 160), fontsize=8, color = new_colr)


        if the_uvb != 'P19' and the_uvb !='FG20':
            txt = 'Q{}'.format(q_val)
            if q_val <= 16:
                ax1.annotate(txt, (n_med+ 0.05, z_med + 0.022), fontsize=8, color = new_colr)
            else:
                ax1.annotate(txt, (n_med - 0.02, z_med + 0.022), fontsize=8, color=new_colr)
        else:
            txt = '{}'.format(the_uvb)
            if the_uvb == 'FG20':
                ax1.annotate(txt, (n_med -0.03, z_med + 0.05), fontsize=8, color = new_colr)
            else:
                ax1.annotate(txt, (n_med -0.03, z_med - 0.1), fontsize=8, color = new_colr)






        print(the_uvb)



    ax1.annotate(r'Hybrid absorber 10$^{5.5}$K' +'\nTrue UVB KS19 (Q18)', xy=(0.11, 0.83), xycoords='axes fraction', fontsize=11)

    ax1.annotate(r'True (log Z, log n$_{\rm H}$) = ' + '({:.0f}, {:.0f})'.format(met, np.log10(den)),
                 xy=(0.11, 0.76), xycoords='axes fraction', fontsize=11)

    x = [-1000, -1000]
    ax1.scatter(x, x, alpha=0.5, s=18, label='using 8 ions', color = 'b', facecolor = 'white')


    lg = ax1.legend(fontsize=11, loc=4)
    lg.get_frame().set_facecolor('none')

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

    add_num0 = 0.42
    ax1.set_xlim(np.log10(den)-add_num0, np.log10(den)+add_num0 + 0.22)
    add_num  = 0.15
    ax1.set_ylim(-1.065, -0.885)
    ax2.set_xlim(np.log10(den)-add_num0, np.log10(den)+add_num0 +0.22)
    ax3.set_xlim(-1.05, -0.93)
    #ax2.set_ylim(0, 174)
    #ax3.set_ylim(0, 174)

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
        make_neon_plots(den=n, met=z, outpath= '/home/vikram/cloudy_run/more_fig/hybrid')
        make_normal_plots(den=n, met=z, outpath= '/home/vikram/cloudy_run/more_fig/hybrid')
