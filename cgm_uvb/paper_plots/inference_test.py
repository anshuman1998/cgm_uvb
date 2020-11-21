import numpy as  np
import corner

path = '/home/vikram/cloudy_run/figures'
file  = path + '/anshuman_try_Q18.npy'

flat_sample = np.load(file)

labels = [r'log n$_{\rm H}$ (cm$^{-3}$)', r'log Z (Z$_{\odot}$)']
truths = [-4, -1]
fig = corner.corner(flat_sample, labels=labels, truths=truths, truth_color = 'b', quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_kwargs={"fontsize": 12})

ndim = 2
axes = np.array(fig.axes).reshape((ndim, ndim))
axes[1,0].annotate ('KS19 (Q18)', xy=(0.03, 0.85), xycoords='axes fraction', fontsize=9)

for ax in [axes[0, 0], axes[1, 0], axes[1,1 ]]:
    ax.tick_params(direction='out', length=5, width=1.1)
    #ax.tick_params(direction='out', which='minor', length=3.5, width=1.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.1)


fig.savefig('inference_Q18.pdf', bbox_inches='tight')

