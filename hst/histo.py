import numpy as  np
import matplotlib as mpl
import astropy.table as tab
import matplotlib.pyplot as plt
import os
import astropy.constants as const
import astropy.units as u
import seaborn as sns

# setting the figure
font = {'family': 'serif', 'weight': 'normal', 'size': 13}
plt.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5

out_fig_name = 'uvb_z01_proposal.pdf'
figure_size = [7, 6]
fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))


d1 = tab.Table.read('pg116_Q14.fits')
d2 = tab.Table.read('pg116_Q20.fits')

sns.distplot(d1['nH'][d1['nH']<-2.75], ax =ax, kde_kws=dict(linewidth=3), label = 'Q14')
sns.distplot(d2['nH'][d2['nH']<-3.25], ax =ax,  kde_kws=dict(linewidth=3), label = 'Q20')

"""
file_names = tab.Table()
for Q in uvb_array:
    q_name = 'Q{}'.format(Q)
    label='KS19 (Q{})'.format(Q)
    if Q ==18:
        ax.plot(12398/uvb['Wave'], uvb[q_name]*4*np.pi*c/(uvb['Wave']*1e-10), label = label, linewidth = 2.5,
                alpha = alpha )
    else:
        ax.plot(12398/uvb['Wave'], uvb[q_name]*4*np.pi*c/(uvb['Wave']*1e-10), label = label, linewidth = 2, linestyle = '--',
                alpha = alpha)

    # 1 angstrom = 12398 eV
"""



ax.legend( loc = 'best', fontsize = 12, ncol=2, handlelength=2.6)
#n_level1 = 'z = {:0.1f} UV Background'.format(z)
ax.annotate ('PG1116 (z = 0.138)', xy=(-4.5, 1.9), fontsize=12)

ax.set_ylabel('PDF')
ax.set_xlabel(r'log n$_{\rm H}$ (cm$^{-3}$)')
#ax.set_xlim (4, 3000)
#ax.set_ylim (4e-8, 8e-5)
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