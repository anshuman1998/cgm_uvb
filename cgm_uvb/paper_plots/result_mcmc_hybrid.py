import numpy as  np
import corner

path = '/home/vikram/cloudy_run/figures/hybrid_precise'

logT_array = [5.0, 5.25, 5.5, 5.75, 6.0, 6.5]

for logT in logT_array:
    uvb = 'KS18'
    Q= 18
    file = path + '/{}_Q{}_logT{:.0f}.npy'.format(uvb, Q, logT*100)
    flat_sample = np.load(file)

    labels = [r'log n$_{\rm H}$ (cm$^{-3}$)', r'log Z (Z$_{\odot}$)']
    truths = [-4, -1]
    fig = corner.corner(flat_sample, labels=labels, quantiles=[0.16, 0.5, 0.84],
                        truths=truths, truth_color='r',
                        show_titles=True, title_kwargs={"fontsize": 12})
    fig.set_size_inches(6.5, 6.5)

    ndim = 2
    axes = np.array(fig.axes).reshape((ndim, ndim))
    uvb_name = uvb
    if uvb == 'KS18':
        uvb_name = 'KS19' + ' (Q{})'.format(Q)
    axes[1, 0].annotate(uvb_name, xy=(0.6, 0.04), xycoords='axes fraction', fontsize=11)
    text_w = 'Warm-hot absorber\n\n'+'at log T (K) = {}\n\n'.format(logT) + r'True (log n$_{\rm H}$, Z) = (-4, -1)' + '\n\nTrue UVB KS19 (Q18)'
    axes[0, 1].annotate(text_w, xy=(0.08, 0.4), xycoords='axes fraction', fontsize=12)


    for ax in [axes[0, 0], axes[1, 0], axes[1, 1]]:
        ax.tick_params(direction='out', length=5, width=1.1)
        # ax.tick_params(direction='out', which='minor', length=3.5, width=1.5)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.1)


    fig.savefig('hybrid_{}_Q{}_logT{:.0f}.pdf'.format(uvb, Q, logT*100), bbox_inches='tight')



