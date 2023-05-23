import numpy as  np
import corner


# for "Photoionized thin model"
path = '/home/vikram/cloudy_run/hst'
file  = path + '/Q17_noshift.npy'

flat_sample = np.load(file)

labels = [r'log n$_{\rm H}$ (cm$^{-3}$)', r'log Z (Z$_{\odot}$)']
truths = [-4, -1]
fig = corner.corner(flat_sample, labels=labels, truths=truths, truth_color = 'r', quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_kwargs={"fontsize": 12})
fig.set_size_inches(6.5, 6.5)

ndim = 2
axes = np.array(fig.axes).reshape((ndim, ndim))
axes[1,0].annotate (r'$\alpha$ = -1.7', xy=(0.6, 0.04), xycoords='axes fraction', fontsize=11)
axes[0,1].annotate ('Photoionized absorber\n'+ 'with ' + r' $\alpha = -1.7$' + '\nUV Background' , xy=(0.08, 0.5), xycoords='axes fraction', fontsize=12)

for ax in [axes[0, 0], axes[1, 0], axes[1,1 ]]:
    ax.tick_params(direction='out', length=5, width=1.1)
    #ax.tick_params(direction='out', which='minor', length=3.5, width=1.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.1)


fig.savefig('inference_Q17_hst.pdf', bbox_inches='tight')
