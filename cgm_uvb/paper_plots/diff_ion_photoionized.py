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

out_fig_name = 'scatter_diff_ion_photo_newfile.pdf'
figure_size = [14, 4]
fig, (ax1, ax2, ax3)  = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.15)


path  = '/home/vikram/cloudy_run/diff_op/photo_NH15'

uvb = ['KS18', 'HM12',  'P19', 'FG20']
uvb_Q = [14, 15, 16, 17, 18, 19, 20]

uvb_models =[]
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

#ks_array = ['14', '15', '16', '17', '18', '19', '20']
#all_uvb = ks_array + ['FG20', 'P19']
#n_H_array = [1e-3, 1e-4, 1e-5]



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

den = 1e-4
met = -1
ion_num = [3, 5,  8]
d = tab.Table.read(path + '/all_combined.fits')

color_list = ['dodgerblue', 'orange', 'green', 'magenta']



for num, col in zip(ion_num, color_list):
    dnum = d[d['n_ions'] == num]
    sort_d = dnum[dnum['true_nH'] == den]
    sort_d = sort_d[sort_d['true_logZ'] == met]

    ax1.scatter(sort_d['n_dmax'], sort_d['z_dmax'], alpha=0.5, label='{} ions'.format(num), s=13)


    binwidth = 0.02

    ax2.hist(sort_d['n_dmax'], bins=np.arange(min(sort_d['n_dmax']), max(sort_d['n_dmax']) + binwidth, binwidth), alpha = 0.75,
             histtype = 'bar', edgecolor= 'k', linewidth = 1.5, density = 1, label = '{} ions'.format(num) )


    ax3.hist(sort_d['z_dmax'], bins=np.arange(min(sort_d['z_dmax']), max(sort_d['z_dmax']) + binwidth, binwidth), alpha = 0.75,
             histtype = 'bar',  edgecolor= 'k', linewidth = 1.5, density = 1, label = '{} ions'.format(num))



    """
    #--------tasted but discarded
    binwidth = 0.025

    p = ax2.hist(n2_all_new, bins=np.arange(min(n2_all_new), max(n2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.35,
                histtype='stepfilled', linewidth=2, density=1, edgecolor ='None',
                label='' ,
                facecolor=col)


    p = ax2.hist(n2_all_new, bins=np.arange(min(n2_all_new), max(n2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.99,
                histtype='step', linewidth=2, density=1, edgecolor =col,
                label='')


    p = ax3.hist(z2_all_new, bins=np.arange(min(z2_all_new), max(z2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.35,
                histtype='stepfilled', linewidth=2, density=1, edgecolor ='None',
                label='' ,
                facecolor=col)

    p = ax3.hist(z2_all_new, bins=np.arange(min(z2_all_new), max(z2_all_new) + binwidth, binwidth)-binwidth/2, alpha=0.99,
                histtype='step', linewidth=2, density=1, edgecolor =col,
                label='' )
    
    """

    #print(np.median(n2_all_new), np.median(z2_all_new), 'ions', num)



ax1.annotate ('Photoionized \n' + 'absorber' , xy=(0.06, 0.65), xycoords='axes fraction', fontsize=12)

ax1.annotate (r'True (log Z, log n$_{\rm H}$) = (-1, -4)' , xy=(0.06, 0.77), xycoords='axes fraction', fontsize=12)

n_med = np.median(sort_d['n_dmax'])
z_med = np.median(sort_d['z_dmax'])

ax2.axvline(n_med, linestyle = '--', color = 'cyan')
ax3.axvline (z_med, linestyle = '--', color = 'cyan')

ax2.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
ax2.annotate (r'$\Delta_{\rm max}$ (log n$_{\rm H}$) =' + '{:.2f}'.format(n_med) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)


ax3.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
ax3.annotate (r'$\Delta_{\rm max}$ (log Z) =' + '{:.2f}'.format(z_med) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)


print(np.median(sort_d['n_dmax']), np.median(sort_d['z_dmax']))

ax1.legend(loc = 'best',  fontsize = 12, ncol=2)

ax1.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
ax1.set_ylabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
#ax1.set_xlim (0.45, 1.02)
#ax1.set_ylim (0.18, 1.02)
#ax2.set_xlim (0.45, 1.05)
#ax3.set_xlim (0.18, 1.05)

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
