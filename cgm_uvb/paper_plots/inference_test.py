import numpy as  np
import corner


# for "Photoionized thin model"
path = '/home/vikram/cloudy_run/figures/2D'
file  = path + '/KS18_Q18.npy'

flat_sample = np.load(file)

labels = [r'log n$_{\rm H}$ (cm$^{-3}$)', r'log Z (Z$_{\odot}$)']
truths = [-4, -1]
fig = corner.corner(flat_sample, labels=labels, truths=truths, truth_color = 'r', quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_kwargs={"fontsize": 12})
fig.set_size_inches(6.5, 6.5)

ndim = 2
axes = np.array(fig.axes).reshape((ndim, ndim))
axes[1,0].annotate ('KS19 (Q18)', xy=(0.6, 0.04), xycoords='axes fraction', fontsize=11)
axes[0,1].annotate ('Photoionized \n'+ 'absorber', xy=(0.08, 0.5), xycoords='axes fraction', fontsize=12)

for ax in [axes[0, 0], axes[1, 0], axes[1,1 ]]:
    ax.tick_params(direction='out', length=5, width=1.1)
    #ax.tick_params(direction='out', which='minor', length=3.5, width=1.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.1)


fig.savefig('inference_Q18.pdf', bbox_inches='tight')



# for "Hybrid model"

path = '/home/vikram/cloudy_run/figures/hybrid'
file  = path + '/KS18_Q18_logT550.npy'

flat_sample = np.load(file)

labels = [r'log n$_{\rm H}$ (cm$^{-3}$)', r'log Z (Z$_{\odot}$)']
truths = [-4, -1]
fig = corner.corner(flat_sample, labels=labels, truths=truths, truth_color = 'r', quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_kwargs={"fontsize": 12})
fig.set_size_inches(6.5, 6.5)


ndim = 2
axes = np.array(fig.axes).reshape((ndim, ndim))
axes[1,0].annotate ('KS19 (Q18)', xy=(0.6, 0.04), xycoords='axes fraction', fontsize=11)
axes[0,1].annotate ('Warm-hot absorber\n' + r'(${\rm T=10^{5.5}}$ K)', xy=(0.08, 0.5), xycoords='axes fraction', fontsize=12)
axes[1,0].set_ylim(-1.35, -0.72)

for ax in [axes[0, 0], axes[1, 0], axes[1,1 ]]:
    ax.tick_params(direction='out', length=5, width=1.1)
    #ax.tick_params(direction='out', which='minor', length=3.5, width=1.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.1)


fig.savefig('inference_Q18_hybrid_T550.pdf', bbox_inches='tight')


