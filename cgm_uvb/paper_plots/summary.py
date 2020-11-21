import numpy as  np
import matplotlib as mpl
import astropy.table as tab
import matplotlib.pyplot as plt
import os

# setting the figure
font = {'family': 'serif', 'weight': 'normal', 'size': 14}
plt.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5

out_fig_name = 'summary.pdf'
figure_size = [12, 5]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(figure_size[0], figure_size[1]))
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0.1)

file_name  =  '/home/vikram/cloudy_run/figures/NH14_out.fits'
d =  tab.Table.read(file_name)

uvb_array= [14, 15, 16, 17, 18, 19, 20]
Q_names = []
for Q in uvb_array:
    q_name = 'Q{}'.format(Q)
    Q_names.append(q_name)

ax1.errorbar( d['nH'], Q_names, xerr = [d['n16'], d['n84']], ls = '', marker ='.', markersize = 14, color = 'b', elinewidth =2, alpha =0.9)
ax2.errorbar( d['Z'], Q_names, xerr = [d['Z16'], d['Z84']], ls = '', marker ='.', markersize = 14, color = 'b', elinewidth=2, alpha = 0.9 )
ax1.axvline(x = -4, linestyle = '--', color = 'k', alpha = 0.7)
ax2.axvline(x = -1, linestyle = '--', color = 'k', alpha = 0.7)

#ax2.set_yticklabels([])


ax1.set_ylabel('UVB models')
ax1.set_xlabel(r'log n$_{\rm H}$ (cm$^{-3}$)')
ax2.set_xlabel(r'log Z (Z$_{\odot}$)')


ax2.set_xlim(-1.25, -0.55)
#ax.set_ylim(4e-8, 5e-4)
#ax.set_xscale('log')
#ax.set_yscale('log')

#ax.legend( loc = 'best', fontsize = 11)
#n_level1 = 'z = {:0.1f}'.format(z)
#ax.annotate (n_level1, xy=(1e4, 2e-7), fontsize=12)
for ax in (ax1, ax2):


    # deco
    ax.tick_params(direction='in', length=6, width=1.5)
    ax.tick_params(direction='in', which='minor', length=3.5, width=1.5)
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    # decorating the plot
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)





fig.savefig(out_fig_name, bbox_inches='tight')