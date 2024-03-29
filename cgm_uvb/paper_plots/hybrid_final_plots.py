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

def make_plot_with_neon_1(den, met, outpath, logT= 5.5, ion_num = [8], uvb_true ='KS18', q_true = 18):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5


    final_path  = outpath
    if not os.path.isdir(final_path):
        os.mkdir(final_path)

    #out_fig_name = outpath + '/hybrid_3plot.pdf'.format(uvb_true, q_true, ion_num[0], met, np.log10(den))
    out_fig_name = outpath +'/hybrid_3plot_logT{:0.0f}.pdf'.format(logT * 100)
    #figure_size = [14, 4]
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))

    figure_size = [5, 10]
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, figsize=(figure_size[0], figure_size[1]))


    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0.2, wspace=0.0)

    path = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'

    #uvb = ['KS18', 'HM12', 'P19', 'FG20']

    uvb = ['KS18',  'P19', 'FG20']

    uvb_Q = [14, 15, 16, 17, 18, 19, 20]

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

    binwidth = 0.006
    binwidth_met = 0.0035

    d = tab.Table.read(path + '/all_combined_logT{:.0f}.fits'.format(logT *100))
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
            ax2.annotate(txt, (n_med+0.014, 90), fontsize=8, color = new_colr)
        else:
            ax2.annotate(txt, (n_med+0.014, 85), fontsize=8, color = new_colr)

        label_name = r' log Z =' + '{:.2f} ->{} Q{}'.format(z_med, the_uvb, q_val)
        ax3.axvline(z_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(z_med)
        if q_val == 16 or q_val == 15 or q_val == 14:
            ax3.annotate(txt, (z_med+0.002, 90), fontsize=8, color = new_colr)
        if the_uvb == 'P19':
            ax3.annotate(txt, (z_med+0.002, 90), fontsize=8, color = new_colr)



        if the_uvb != 'P19' and the_uvb !='FG20':
            txt = 'Q{}'.format(q_val)
            if q_val <= 16:
                ax1.annotate(txt, (n_med -0.02, z_med + 0.022), fontsize=8, color = new_colr)
            else:
                ax1.annotate(txt, (n_med - 0.02, z_med + 0.022), fontsize=8, color=new_colr)
        else:
            txt = '{}'.format(the_uvb)
            if the_uvb == 'FG20':
                ax1.annotate(txt, (n_med -0.03, z_med + 0.04), fontsize=8, color = new_colr)
            else:
                ax1.annotate(txt, (n_med -0.025, z_med + 0.03), fontsize=8, color = new_colr)

        print(the_uvb)



    ax1.annotate('Warm-hot absorber log T (K) ={:0.1f}'.format(logT) +'\nTrue UVB KS19 (Q18)', xy=(0.11, 0.83), xycoords='axes fraction', fontsize=11)

    ax1.annotate(r'True (log Z, log n$_{\rm H}$) = ' + '({:.0f}, {:.0f})'.format(met, np.log10(den)),
                 xy=(0.11, 0.76), xycoords='axes fraction', fontsize=11)

    x = [-1000, -1000]
    ax1.scatter(x, x, alpha=0.5, s=18, label='Using 8 ions (including Ne VIII)', color = 'b', facecolor = 'white')

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

    add_num0 = 0.31
    ax1.set_xlim(np.log10(den)-add_num0, np.log10(den)+add_num0 + 0.25)
    add_num  = 0.15
    ax1.set_ylim(-1.12, -0.88)
    ax2.set_xlim(np.log10(den)-add_num0, np.log10(den)+add_num0 +0.25)
    ax3.set_xlim(-1.085, -0.93)
    #ax2.set_ylim(0, 174)
    #ax3.set_ylim(0, 174)
    ax2.set_ylim(0, 99)
    ax3.set_ylim(0, 99)

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

#--------------------------for n = 1e-3
def make_plot_with_neon_2(den, met, outpath, logT= 5.5, ion_num = [8], uvb_true ='KS18', q_true = 18):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5


    final_path  = outpath
    if not os.path.isdir(final_path):
        os.mkdir(final_path)

    #out_fig_name = outpath + '/hybrid_4plot.pdf'.format(uvb_true, q_true, ion_num[0], met, np.log10(den))
    out_fig_name = outpath + '/hybrid_4plot_logT{:0.0f}.pdf'.format(logT*100)

    #figure_size = [14, 4]
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))

    figure_size = [5, 10]
    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, figsize=(figure_size[0], figure_size[1]))


    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0.2, wspace=0.0)

    path = '/home/vikram/cloudy_run/diff_op/hybrid_NH15'

    #uvb = ['KS18', 'HM12', 'P19', 'FG20']

    uvb = ['KS18',  'P19', 'FG20']

    uvb_Q = [14, 15, 16, 17, 18, 19, 20]

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

    binwidth = 0.003
    binwidth_met = 0.001

    d = tab.Table.read(path + '/all_combined_logT{:.0f}.fits'.format(logT *100))
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

    dd = 0
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
        if the_uvb == 'P19' or the_uvb == 'FG20':
            ax2.annotate(txt, (n_med + 0.007, 85), fontsize=8, color=new_colr)
        else:
            if q_val == 19:
                ax2.annotate(txt, (n_med + 0.007, 80), fontsize=8, color=new_colr)
            else:
                ax2.annotate(txt, (n_med + 0.007, 90), fontsize=8, color=new_colr)

        label_name = r' log Z =' + '{:.2f} ->{} Q{}'.format(z_med, the_uvb, q_val)
        ax3.axvline(z_med, label=label_name, linewidth = 0.8, dashes=[5, 7], color = new_colr, zorder =50)
        txt = '{:.2f}'.format(z_med)
        if the_uvb =='FG20' or q_val ==14:
            ax3.annotate(txt, (z_med + 0.0002, 85), fontsize=8, color=new_colr)
        if the_uvb == 'KS18' and q_val ==18:
            ax3.annotate(txt, (z_med + 0.0002, 85), fontsize=8, color=new_colr)


        #dd = dd+8




        if the_uvb != 'P19' and the_uvb !='FG20':
            txt = 'Q{}'.format(q_val)
            if q_val <= 17 or q_val ==19:
                ax1.annotate(txt, (n_med -0.01, z_med + 0.007), fontsize=8, color = new_colr)
            if q_val == 20:
                ax1.annotate(txt, (n_med -0.01, z_med + 0.005), fontsize=8, color = new_colr)
            if q_val == 18:
                ax1.annotate(txt, (n_med - 0.01, z_med - 0.004), fontsize=8, color=new_colr)
        else:
            txt = '{}'.format(the_uvb)
            if the_uvb == 'FG20':
                ax1.annotate(txt, (n_med -0.01, z_med - 0.006), fontsize=8, color = new_colr)
            else:
                ax1.annotate(txt, (n_med -0.01, z_med + 0.004), fontsize=8, color = new_colr)


        print(the_uvb)




    ax1.annotate('Warm-hot absorber log T (K) = {:.1f}'.format(logT) +'\nTrue UVB KS19 (Q18)', xy=(0.11, 0.83), xycoords='axes fraction', fontsize=11)

    ax1.annotate(r'True (log Z, log n$_{\rm H}$) = ' + '({:.0f}, {:.0f})'.format(met, np.log10(den)),
                 xy=(0.11, 0.76), xycoords='axes fraction', fontsize=11)

    x = [-1000, -1000]
    ax1.scatter(x, x, alpha=0.5, s=18, label='Using 8 ions (including Ne VIII)', color = 'b', facecolor = 'white')

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

    add_num0 = 0.1
    ax1.set_xlim(-3.13, -2.57)
    add_num  = 0.15
    ax1.set_ylim(-1.03, -0.96)
    ax2.set_xlim(-3.13, -2.57)
    ax3.set_xlim(-1.013, -0.983)
    ax2.set_ylim(0, 99)
    ax3.set_ylim(0, 99)

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


Zarray = [-1]
Narray = [1e-5]
for n in Narray:
    for z in Zarray:
        make_plot_with_neon_1(den=n, met=z, logT=5.75, outpath= '/home/vikram/cloudy_run/diff_op/fig')


Zarray = [-1]
Narray = [1e-3]
for n in Narray:
    for z in Zarray:
        make_plot_with_neon_2(den=n, met=z, logT= 5.75, outpath= '/home/vikram/cloudy_run/diff_op/fig')