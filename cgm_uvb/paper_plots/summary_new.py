import numpy as  np
import matplotlib as mpl
import astropy.table as tab
import matplotlib.pyplot as plt
import os

def make_summary_plot_new(file_name_fits, lsf_file_name, fig_name, add_fits = None, add_lsf = None,  add_log_lsf = False,
                          add_label = '',  x1lim= [-4.3, -3.45], x2lim = [-1.3, -0.65], fig_text = '', loglsf = False):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5

    figure_size = [12, 5]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(figure_size[0], figure_size[1]))
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.05)

    d = tab.Table.read(file_name_fits)
    dlsf = tab.Table.read(lsf_file_name)

    lab = 'MCMC'
    num = np.arange(len(d['uvb']))
    ax1.errorbar(d['nH'][:-1], num[:-1], xerr=[d['n16'][:-1], d['n84'][:-1]], ls='', marker='.', markersize=16,
        color='b', elinewidth=2, alpha=0.8, label=lab)
    if add_fits != None:
        dn = tab.Table.read(add_fits)
        lab = add_label
        ax1.errorbar(dn['nH'][:-1], num[:-1]+0.1, xerr=[dn['n16'][:-1], dn['n84'][:-1]], ls='', marker='.', markersize=16,
                     color='gold', elinewidth=2, alpha=0.99, label=lab, zorder  =-1)
        for i in num:
            print('Delta N= ',  d['nH'][i]- dn['nH'][i], 'for', d['uvb'][i],  d['Q'][i] )


    lab = 'Least Squares'
    if loglsf:
        ax1.scatter(np.log10(dlsf['nH'][:-1]), num[:-1], marker = 'x', color ='r', linewidth = 2, s = 300,  alpha =0.7, label=lab)
    else:
        ax1.scatter(dlsf['nH'][:-1], num[:-1], marker='x', color='r', linewidth=2, s=300, alpha=0.7, label=lab)

    lab = 'MCMC'
    ax2.errorbar(d['Z'][:-1], num[:-1], xerr=[d['Z16'][:-1], d['Z84'][:-1]], ls='', marker='.', markersize=16,
        color='b', elinewidth=2, alpha=0.7)
    if add_fits:
        dn = tab.Table.read(add_fits)
        ax2.errorbar(dn['Z'][:-1], num[:-1] + 0.1, xerr=[dn['Z16'][:-1], dn['Z84'][:-1]], ls='', marker='.', markersize=16,
                     color='gold', elinewidth=2, alpha=0.99, label=add_label, zorder=-1)

    lab = 'Least Squares'
    if loglsf:
        ax2.scatter(np.log10(dlsf['Z'][:-1]), num[:-1], marker = 'x', color ='r', linewidth = 2, s = 300,  alpha =0.7)
    else:
        ax2.scatter(dlsf['Z'][:-1], num[:-1], marker='x', color='r', linewidth=2, s=300, alpha=0.7)

    ax1.axvline(x=-4, linestyle='--', color='k', alpha=0.7)
    ax2.axvline(x=-1, linestyle='--', color='k', alpha=0.7)

    # ax1.tick_params(axis='x', which='minor', bottom=False)
    ax1.set_yticks(num[:-1])

    uvb_column = ['Q14', 'Q15', 'Q16', 'Q17', 'Q18', 'Q19', 'Q20', 'P19', 'FG20', 'HM12']
    ax1.set_yticklabels(uvb_column[:-1])
    ax2.set_yticklabels([])

    ax1.set_ylabel('UVB models')
    ax1.set_xlabel(r'log n$_{\rm H}$ (cm$^{-3}$)')
    ax2.set_xlabel(r'log Z (Z$_{\odot}$)')

    ax1.set_xlim(x1lim[0], x1lim[1])
    ax2.set_xlim(x2lim[0], x2lim[1])
    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax2.annotate(fig_text, xy=(0.65, 0.8), xycoords='axes fraction', fontsize=14)

    # change the marker size manually for both lines
    # lgnd.legendHandles[0]._legmarker.set_markersize(6)
    # lgnd.legendHandles[1]._legmarker.set_markersize(100)
    # n_level1 = 'z = {:0.1f}'.format(z)
    # ax.annotate (n_level1, xy=(1e4, 2e-7), fontsize=12)
    for ax in (ax1, ax2):

        if ax ==ax1:
            ax.legend(loc='best', fontsize=14)

        # deco
        ax.tick_params(direction='in', length=6, width=1.5)
        ax.tick_params(direction='in', which='minor', length=3.5, width=1.5)
        # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        # decorating the plot
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.5)

    fig.savefig(fig_name, bbox_inches='tight')

    return



lsf_file_name = '/home/vikram/cloudy_run/figures/2D/NH15_log_lsf_out.fits'
file_name_phot = '/home/vikram/cloudy_run/figures/2D/NH15_metal_2D.fits'
out_fig_name = 'summary_photoionized.pdf'
txt = 'Photoionized\nabsorber'
make_summary_plot_new(file_name_fits= file_name_phot, lsf_file_name =lsf_file_name, fig_name= out_fig_name,
                      fig_text= txt)

file_name  =  '/home/vikram/cloudy_run/figures/hybrid/NH15_hybrid_logT550.fits'
lsf_file_name = '/home/vikram/cloudy_run/figures/hybrid/NH15_log_lsf_hybrid_T550.fits'
out_fig_name = 'summary_hybrid.pdf'
txt = 'Warm-hot absorber\n' + r'(T = 10$^{5.5}$ K)'
make_summary_plot_new(file_name_fits= file_name, lsf_file_name =lsf_file_name, add_fits= file_name_phot,
                      fig_name= out_fig_name, loglsf= True, add_label= 'photoionized', fig_text= txt)


"""
#file_name  =  '/home/vikram/cloudy_run/figures/hybrid/NH15_hybrid_logT500.fits'
#lsf_file_name = '/home/vikram/cloudy_run/figures/hybrid/NH15_log_lsf_hybrid_T500.fits'
#out_fig_name = 'summary_hybrid.pdf'
#make_summary_plot(file_name_fits= file_name, lsf_file_name =lsf_file_name, fig_name= out_fig_name, loglsf= True)

#---------------


file_name  =  '/home/vikram/cloudy_run/figures/rescaled/rescaled_NH15_metal_2D.fits'
lsf_file_name = '/home/vikram/cloudy_run/figures/rescaled/rescaled_NH15_log_lsf_out.fits'
out_fig_name = 'summary_rescaled.pdf'
make_summary_plot(file_name_fits= file_name, lsf_file_name =lsf_file_name, fig_name= out_fig_name, loglsf= False)


"""

print('Hybrid')
file_name  =  '/home/vikram/cloudy_run/figures/rescaled_hybrid/NH15_hybrid_logT550.fits'
lsf_file_name = '/home/vikram/cloudy_run/figures/rescaled_hybrid/NH15_log_lsf_hybrid_T550.fits'
add_fits =  '/home/vikram/cloudy_run/figures/hybrid/NH15_hybrid_logT550.fits'
add_lsf = '/home/vikram/cloudy_run/figures/hybrid/NH15_log_lsf_hybrid_T550.fits'
out_fig_name = 'summary_rescaled_hybrid_new.pdf'
txt = 'Warm-hot absorber\n' + r'(T = 10$^{5.5}$ K)' +'\nrescaled UVB'
make_summary_plot_new(file_name_fits= file_name, lsf_file_name =lsf_file_name, fig_name= out_fig_name,
                      add_lsf= add_lsf, loglsf= True, fig_text= txt,
                      add_fits  =  add_fits, add_label= 'original')

print('Photoionized')
file_name  =  '/home/vikram/cloudy_run/figures/rescaled/rescaled_NH15_metal_2D.fits'
lsf_file_name = '/home/vikram/cloudy_run/figures/rescaled/rescaled_NH15_log_lsf_out.fits'
add_fits =  '/home/vikram/cloudy_run/figures/2D/NH15_metal_2D.fits'
add_lsf = '/home/vikram/cloudy_run/figures/2D/NH15_log_lsf_out.fits'
out_fig_name = 'summary_rescaled_photoinized_new.pdf'
txt = 'Photoionized\nabsorber\nrescaled UVB'
make_summary_plot_new(file_name_fits= file_name, lsf_file_name =lsf_file_name, fig_name= out_fig_name,
                      add_lsf= add_lsf, loglsf= False, fig_text= txt,
                      add_fits  =  add_fits, add_label= 'original')