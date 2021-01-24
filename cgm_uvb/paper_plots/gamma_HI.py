import astropy.table as t
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
from matplotlib.ticker import ScalarFormatter
from astropy.cosmology import Planck15 as planck
import scipy.optimize as sciopt


# setting the figure
font = {'family': 'serif', 'weight': 'normal', 'size': 14}
plt.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.5

out_fig_name = 'uvb_z02.pdf'
figure_size = [7, 6]
fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))

alpha=0.9
figname='Gamma_HI.pdf'

filename='/home/vikram/Work/centOS/UVB/KS2016/uvb/models/data/plots/gama/becker2013.txt'
data=np.loadtxt(filename)
z=data[:,0]
g=data[:,1]
g2=data[:,2]
g1=data[:,3]
#ax.errorbar(z, 10.0**g*1e-12, yerr=[(10.0**(g)-10.0**(g+g1))*1e-12, (10.0**(g+g2)-10.0**(g))*1e-12,], marker='s', label='Becker et al. 2013', markersize=8, ls='', c='darkorchid', alpha=0.9, elinewidth=2)

filename='/home/vikram/Work/centOS/UVB/KS2016/uvb/models/data/plots/gama/bolton07.txt'
data=np.loadtxt(filename)
z=data[:,0]
g=data[:,1]
g2=data[:,2]
g1=data[:,3]
#ax.errorbar(z, g*1e-12, yerr=[-1.0*g1*1e-12, g2*1e-12,], marker='v', label='Bolton et al. 2007', markersize=11, ls='', c='green', alpha=0.9, elinewidth=2)

x=0.1
y=0.175e-12
l='Kollmeier et al. 2014'
#ax.scatter(x, y, label=l, marker='*', s=300, c='blue')

uvb_array= [14, 15, 16, 17, 18, 19, 20]

al = 0.9
for Q in uvb_array:
  q_name = 'Q{}'.format(Q)
  label = 'KS19 ({})'.format(q_name)
  uvb = np.loadtxt('/home/vikram/Downloads/KS_2018_EBL/parameters_Q{}.txt'.format(Q))
  z = uvb[:, 0]
  g = uvb[:, 1]
  #ax.plot(z, g, color='k', label=label, linestyle='-', dashes=(5, 4), linewidth=2, alpha= al)
  ax.plot(z, g, label=label, linewidth=2, alpha= al)

  gamma_file = '/home/vikram/cgm_uvb/cgm_uvb/paper_plots/gamma_HI_files/' + 'Gamma_HI_KS18_Q{}.fits'.format(Q)
  data = t.Table.read(gamma_file)
  ax.plot(data['z'], data['g'], linestyle='-', dashes=(5, 4), linewidth=2, alpha= al, color = 'k')



#ax.plot(z, g, color='k', label=r'Khaire & Srianand 2018 ', linewidth=2, alpha=0.8)

uvb=np.loadtxt('/home/vikram/Work/data_literature/gamma_h1/HM_gama.dat')
z=uvb[:,0]
g=uvb[:,1]
ax.plot(z, g, color='magenta', label='Haardt & Madau 2012', linestyle='-.', linewidth=2, dashes=(5, 4, 1, 2), alpha=0.8)
#ax.plot(z, g, color='magenta', label='Haardt & Madau 2012', alpha=0.8)

filename='/home/vikram/Work/ucsb/chi_square/diagonal_best_fit_gamma.txt'
data=np.loadtxt(filename)
z=data[:,0]
g1=data[:,1]
g=data[:,2]
g2=data[:,3]
i=2
#ax.errorbar(z, g, (g2-g), marker='.', markersize=13, ls='', c='red', alpha=0.5, elinewidth=2)
#print(z, g[i])
ax.errorbar(z, g, yerr=[5.75*(g-g1), 5.75*(g2-g),], marker='.', label='Khaire et al. 2019', markersize=20, ls='',
            c='b', alpha=0.9, elinewidth=2.9, zorder=11, capsize=4, capthick=2.2)

# Gaikwad et al 2018
pz=np.array([0.11, 0.21, 0.31, 0.41])
pg=np.array([0.066, 0.1, 0.145, 0.210])*1e-12
pe=np.array([0.015, 0.021, 0.037, 0.052])*1e-12
t = ['g', 'orange', 'b', 'm', 'r']
l='Gaikwad et al. 2017'
ax.errorbar(pz, pg, pe, marker='D', label=l, markersize=8, ls='', c='k',  alpha=0.8, elinewidth=2.5, zorder = 10, capsize=4, capthick=2.2)


l='Fumagalli + in prep'
z=[0.004,]
g=[6.39e-14,]
err1=[2.73e-14,]
err2=[2.91e-14,]
ax.errorbar(z, g, yerr=[err1, err2,],  marker='s', label=l, markersize=8, zorder=11,  c='green',  alpha=alpha,  elinewidth=2.5, clip_on=False)

ax.set_yscale('log')
ax.set_ylabel(r'$\Gamma_{\rm H \, I}$ ( s$^{-1} )$')
ax.set_xlabel('Redshift')
ax.legend( loc = 'best', fontsize = 11, handlelength=2.8)


ax.set_xlim(-0.005, 0.55)
ax.set_ylim(1e-14, 4e-13)
ax.tick_params(direction='in', length=7, width=1.7)
ax.tick_params(direction='in', which='minor', length=4, width=1.7)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
#ax.ticklabel_format(axis='y', style='sci')
#ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))



for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.7)

ax.legend(loc='lower right')

fig.tight_layout(rect=[-0.03, -0.03, 1.02, 1.02])

fig.savefig(figname, bbox_inches='tight')

plt.show()

