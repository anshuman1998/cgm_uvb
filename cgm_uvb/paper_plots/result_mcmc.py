import numpy as  np
import corner

path = '/home/vikram/cloudy_run/figures/2D'

uvb_array = ['KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'KS18', 'P19', 'FG20', 'HM12']
Q_array= [14, 15, 16, 17, 18, 19, 20, 18, 18, 18]

for Q, uvb in zip(Q_array, uvb_array):
    file = path + '/{}_Q{}.npy'.format(uvb, Q)
    flat_sample = np.load(file)

    labels = [r'log n$_{\rm H}$ (cm$^{-3}$)', r'log Z (Z$_{\odot}$)']
    fig = corner.corner(flat_sample, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})

    ndim = 2
    axes = np.array(fig.axes).reshape((ndim, ndim))
    if uvb == 'KS18':
        legend_name = 'KS19 (Q{})'.format(Q)
    else:
        legend_name = '{}'.format(uvb)
    axes[1, 0].annotate(legend_name, xy=(0.03, 0.85), xycoords='axes fraction', fontsize=9)

    for ax in [axes[0, 0], axes[1, 0], axes[1, 1]]:
        ax.tick_params(direction='out', length=5, width=1.1)
        # ax.tick_params(direction='out', which='minor', length=3.5, width=1.5)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(1.1)

    fig.savefig('result_{}_Q{}.pdf'.format(uvb, Q), bbox_inches='tight')



