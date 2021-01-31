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

out_fig_name = 'scatter_hybrid.pdf'
figure_size = [7, 6]
fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))

path  = '/home/vikram/cgm_uvb/cgm_uvb/paper_plots/more_models'

#d  = tab.Table.read(path + '/diff_res_2ions_new.txt', format = 'ascii')
ks_array = ['14', '15', '16', '17', '18', '19', '20']
all_uvb = ks_array + ['FG20', 'P19']
n2_ks = []
z2_ks = []
n2_all = []
z2_all = []

ion_num = [3, 4, 5, 6]

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

for num in ion_num:
    d = tab.Table.read(path + '/diff_res_{}ions_newhyb.txt'.format(num), format='ascii')
    print(num,  len(d))
    n2_all_new = []
    z2_all_new = []
    for m, n, Z, in zip(d['Q'], d['Max_diff_nH'], d['Max_diff_Z']):
        for uvb in all_uvb:
            if uvb == m:
                n2_all_new.append(n)
                z2_all_new.append(Z)
                n2_all.append(n)
                z2_all.append(Z)
    print(len(n2_all))
    ax.scatter(n2_all_new, z2_all_new, alpha = 0.7, label = '{} ions'.format(num), s = 3)

print(np.median(n2_all), np.median(z2_all))

ax.legend(loc = 'best')

ax.set_ylabel(r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)')
ax.set_xlabel(r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)')
#ax.set_xlim (0, 1.5)
#ax.set_ylim (0, 1.5)
#ax.set_xscale('log')
#ax.set_yscale('log')

#deco
ax.tick_params(direction='in', length=7, width=1.7)
ax.tick_params(direction='in', which='minor', length=4, width=1.7)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
# decorating the plot
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.7)
    ax.spines[axis].set_color('k')

fig.savefig(out_fig_name, bbox_inches='tight')


import corner
import numpy as np

labels = [r'$\Delta_{\rm max}$ log n$_{\rm H}$ (cm $^{-3}$)', r'$\Delta_{\rm max}$ log Z(Z$_{\odot}$)']

samples = (np.vstack([n2_all, z2_all])).T

figure = corner.corner(samples,  labels=labels, quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig('test_s.pdf', bbox_inches='tight')


