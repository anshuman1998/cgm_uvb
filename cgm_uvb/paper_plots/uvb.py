import numpy as  np
import matplotlib as mpl
import astropy.table as tab
import matplotlib.pyplot as plt
import os
import astropy.constants as const
import astropy.units as u

# setting the figure
font = {'family': 'serif', 'weight': 'normal', 'size': 13}
plt.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5

out_fig_name = 'uvb_z02.pdf'
figure_size = [7, 6]
fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))
c = (const.c).to(u.m/u.s).value # speed of light in m/s
z = 0.2
path = os.getcwd() + '/ks19_ebl'
file_name  =  path + '/KS19_EBL_z_{:.1f}.fits'.format(z)
uvb =  tab.Table.read(file_name)

uvb_array= [14, 15, 16, 17, 18, 19, 20]

file_names = tab.Table()
for Q in uvb_array:
    q_name = 'Q{}'.format(Q)
    label= 'KS19 ({})'.format(q_name)
    ax.plot(12398/uvb['Wave'], uvb[q_name]*4*np.pi*c/(uvb['Wave']*1e-10), label = label, linewidth = 2)
    # 1 angstrom = 12398 eV

fg20_path = os.getcwd() + '/fg20_fits_files'
file_name  =  fg20_path + '/FG20_EBL_z_{:.2f}.fits'.format(z)
uvb =  tab.Table.read(file_name)
label= 'FG20'
ax.plot(12398/uvb['Wave'], uvb['Jnu'] * 4 * np.pi * c / (uvb['Wave'] * 1e-10), label=label, linewidth=2, linestyle  = '--')
# 1 ryd = 13.6057 eV

ax.legend( loc = 'best', fontsize = 11)
n_level1 = 'z = {:0.1f}'.format(z)
ax.annotate (n_level1, xy=(1e4, 2e-7), fontsize=12)

ax.set_ylabel(r'4$\pi$$\nu$J$_{\nu}$ (ergs s$^{-1}$ cm $^{-2}$)')
ax.set_xlabel(r'Energy (eV)')
ax.set_xlim (1, 100000)
ax.set_ylim (4e-8, 5e-4)
ax.set_xscale('log')
ax.set_yscale('log')

#deco
ax.tick_params(direction='in', length=6, width=1.5)
ax.tick_params(direction='in', which='minor', length=3.5, width=1.5)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
# decorating the plot
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

fig.savefig(out_fig_name, bbox_inches='tight')