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

out_fig_name = 'scatter_diff_ion_hybrid_NH4_bimodal_with_Ne8.pdf'
figure_size = [14, 4]
fig, (ax1, ax2, ax3)  = plt.subplots(1, 3, figsize=(figure_size[0], figure_size[1]))
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.15)


path  = '/home/vikram/cgm_uvb/cgm_uvb/paper_plots/more_models'

#d  = tab.Table.read(path + '/diff_res_2ions_new.txt', format = 'ascii')
ks_array = ['14', '15', '16', '17', '18', '19', '20']
all_uvb = ks_array + ['FG20', 'P19']
all_uvb = [18]
n_H_array = [1e-3, 1e-4, 1e-5]
n2_ks = []
z2_ks = []
nall_low = []
zall_low = []
nall_high = []
zall_high = []

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

den = 1e-5
ion_num = [9]
color_list = ['dodgerblue']
#, 'grey', 'green', 'magenta',  'cyan', 'gold']

left_hand = False
num = ion_num[0]
d = tab.Table.read(path + '/full_{}ions_t550.txt'.format(num), format='ascii')
ions_used_coloum  = d['Ions_Used']
flag = []
for i in range(len(d)):
    ions = ions_used_coloum[i].split('_')
    ne8 = False
    for j in ions:
        if j=='Ne+7':
            ne8= True
    if ne8 == True:
        flag.append(1)
    else:
        flag.append(0)

flag = np.array(flag)
for num, colr, in zip(ion_num, color_list):
    d = tab.Table.read(path + '/full_{}ions_t550.txt'.format(num), format='ascii')

    dlow = d [flag==1]
    dhigh = d [flag==0]
    print(len(dlow), len(dhigh))


    n2_all_new = []
    z2_all_new = []
    n2_all_new_high = []
    z2_all_new_high = []
    for m, n, Z, in zip(dlow['Q'][dlow['true_nH'] == den], dlow['Max_diff_nH'][dlow['true_nH'] == den],
                        dlow['Max_diff_Z'][dlow['true_nH'] == den]):
        for uvb in all_uvb:
            if uvb == m:
                n2_all_new.append(n)
                z2_all_new.append(Z)
                nall_low.append(n)
                zall_low.append(Z)

    for m, n, Z, in zip(dhigh['Q'][dhigh['true_nH'] == den], dhigh['Max_diff_nH'][dhigh['true_nH'] == den],
                        dhigh['Max_diff_Z'][dhigh['true_nH'] == den]):
        for uvb in all_uvb:
            if uvb == m:
                n2_all_new_high.append(n)
                z2_all_new_high.append(Z)
                nall_high.append(n)
                zall_high.append(Z)
    #print(len(n2_all))

    #------
    binwidth = 0.02

    ax1.scatter(n2_all_new, z2_all_new, alpha=0.5, label='{} ions'.format(num), s=13, color = colr)

    ax2.hist(n2_all_new, bins=np.arange(np.min(n2_all_new), np.max(n2_all_new) + binwidth, binwidth), alpha = 0.75,
             histtype = 'bar', edgecolor= 'k', linewidth = 1.5, density = 1, label = '{} ions'.format(num), color = colr )
    ax3.hist(z2_all_new, bins=np.arange(np.min(z2_all_new), np.max(z2_all_new) + binwidth, binwidth), alpha = 0.8,
             histtype = 'bar',  edgecolor= 'k', linewidth = 1.5, density = 1, label = '{} ions'.format(num), color  = colr)



    ax1.scatter(n2_all_new_high, z2_all_new_high, alpha=0.5, label='', s=13, color = 'r')

    ax2.hist(n2_all_new_high, bins=np.arange(np.min(n2_all_new_high), np.max(n2_all_new_high) + binwidth, binwidth), alpha = 0.75,
             histtype = 'bar', edgecolor= 'k', linewidth = 1.5, density = 1, label = '', color = 'r')
    ax3.hist(z2_all_new_high, bins=np.arange(np.min(z2_all_new_high), np.max(z2_all_new_high) + binwidth, binwidth), alpha = 0.8,
             histtype = 'bar',  edgecolor= 'k', linewidth = 1.5, density = 1, label = '', color  = 'r')


    print(np.median(n2_all_new), np.median(z2_all_new), 'for ions', num)

ax1.annotate ('Hybrid \n' + r'absorber (10$^{5.5}$ K) ' , xy=(0.06, 0.75), xycoords='axes fraction', fontsize=12)
ax1.annotate (r'True (log Z, log n$_{\rm H}$) = (-1, -4)' , xy=(0.06, 0.67), xycoords='axes fraction', fontsize=12)


ax2.axvline(np.median(nall_low), linestyle = '--', color = 'blue', zorder = 20)
ax3.axvline (np.median(zall_low), linestyle = '--', color = 'blue', zorder = 20)
ax2.axvline(np.median(nall_high), linestyle = '--', color = 'red', zorder = 20)
ax3.axvline (np.median(zall_high), linestyle = '--', color = 'red', zorder = 20)

ax2.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
ax2.annotate (r'$\Delta_{\rm max}$ (log n$_{\rm H}$)' + '= {:.2f}'.format(np.median(nall_low)) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)
ax2.annotate (r'$\Delta_{\rm max}$ (log n$_{\rm H}$)' + '= {:.2f}'.format(np.median(nall_high)) , xy=(0.6, 0.73), xycoords='axes fraction', fontsize=11)


ax3.annotate (r'Median' , xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
ax3.annotate (r'$\Delta_{\rm max}$' +' (log Z) = {:.2f}'.format(np.median(zall_low)) , xy=(0.6, 0.83), xycoords='axes fraction', fontsize=11)
ax3.annotate (r'$\Delta_{\rm max}$' +' (log Z) = {:.2f}'.format(np.median(zall_high)) , xy=(0.6, 0.73), xycoords='axes fraction', fontsize=11)


##print(np.median(n2_all), np.median(z2_all))

ax1.legend(loc = 'best',  fontsize = 12, ncol=2)

ax1.set_xlabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
ax1.set_ylabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
#ax1.set_xlim (0.5, 1.4)
#ax1.set_ylim (0.25, 1.02)
ax2.set_xlim (0.45, 1.07)
ax3.set_xlim (-0.05, 0.55)

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