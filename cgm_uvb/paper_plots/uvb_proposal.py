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

out_fig_name = 'uvb_z01_proposal.pdf'
figure_size = [7, 6]
fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))
c = (const.c).to(u.m/u.s).value # speed of light in m/s
z = 0.1
alpha = 0.99
path = os.getcwd() + '/ks19_ebl'
file_name  =  path + '/KS19_EBL_z_{:.1f}.fits'.format(z)
uvb =  tab.Table.read(file_name)

uvb_array= [14, 15, 16, 17, 18, 19, 20]

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

fg20_path = os.getcwd() + '/fg20_fits_files'
file_name  =  fg20_path + '/FG20_EBL_z_{:.2f}.fits'.format(z)
uvb =  tab.Table.read(file_name)
label= 'FG20'
ax.plot(12398/uvb['Wave'], uvb['Jnu'] * 4 * np.pi * c / (uvb['Wave'] * 1e-10), label=label, linewidth=2,
        linestyle  = '-.',alpha = alpha, c= 'grey')

#fg20_path = os.getcwd() + '/p19_ebl'
#file_name  =  fg20_path + '/P19_EBL_z_{:.2f}.fits'.format(z)
#uvb =  tab.Table.read(file_name)
#label= 'P19'
#ax.plot(12398/uvb['Wave'], uvb['Jnu'] * 4 * np.pi * c / (uvb['Wave'] * 1e-10), label=label, linewidth=3,
#        linestyle  = '-.', alpha = alpha, c= 'cyan')

# 1 ryd = 13.6057 eV


#for ionization potentials
ax.vlines(x =  23.383143, ymin = 1.5e-06, ymax = 2.5e-06,colors = 'grey',linestyles ="solid", lw= 1.0)
ax.text(21.08,2.8e-06 , "C II", fontsize=8)
ax.vlines(x =  64.49352, ymin = 1.5e-06, ymax = 2.5e-06,colors = 'grey',linestyles ="solid",lw= 1.0)
ax.text(59.2,2.8e-06 , "C IV", fontsize=8)
ax.vlines(x =  29.60125, ymin = 1.5e-06, ymax = 3.0e-06,colors = 'grey',linestyles ="solid",  lw= 1.0)
ax.text(26.20125,3.3e-06 , "N II", fontsize=8)
ax.vlines(x =  97.8901, ymin = 1.5e-06, ymax = 2.5e-06,colors = 'grey',linestyles ="solid",lw= 1.0)
ax.text(91.1901,2.8e-06 , "N V", fontsize=8)
ax.vlines(x =  138.1189,ymin = 1.5e-06, ymax = 3.0e-06,colors = 'grey',linestyles ="solid", lw= 1.0)
ax.text(126.7189,3.3e-06 , "O VI", fontsize=8)
ax.vlines(x =  16.34585, ymin = 1.5e-06, ymax = 3.0e-06,colors = 'grey',linestyles ="solid",lw= 1.0)
ax.text(14.84585,3.3e-06 , "Si II", fontsize=8)
ax.vlines(x =   33.49300, ymin = 1.5e-06, ymax = 2.5e-06,colors = 'grey',linestyles ="solid",lw= 1.0)
ax.text(30.49300,2.8e-06 , "Si III", fontsize=8)
ax.vlines(x =   45.14179, ymin = 1.5e-06, ymax = 3.0e-06,colors = 'grey',linestyles ="solid",lw= 1.0)
ax.text(41.14179,3.3e-06 , "Si IV", fontsize=8)
ax.vlines(x =   239.0970, ymin = 1.5e-06, ymax = 2.5e-06,colors = 'grey',linestyles ="solid",lw= 1.0)
ax.text(201.0970,2.8e-06 , "Ne VIII", fontsize=8)
ax.vlines(x =   367.489, ymin = 1.5e-06, ymax = 3.0e-06,colors = 'grey',linestyles ="solid",lw= 1.0)
ax.text(332.489,3.3e-06 , "Mg X", fontsize=8)

ax.legend( loc = 'best', fontsize = 12, ncol=2, handlelength=2.6)
n_level1 = 'z = {:0.1f} UV Background'.format(z)
ax.annotate (n_level1, xy=(0.7e2, 0.8e-5), fontsize=12)

ax.set_ylabel(r'4$\pi$$\nu$J$_{\nu}$ (ergs s$^{-1}$ cm $^{-2}$)')
ax.set_xlabel(r'Energy (eV)')
ax.set_xlim (4, 3000)
ax.set_ylim (4e-8, 8e-5)
ax.set_xscale('log')
ax.set_yscale('log')

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