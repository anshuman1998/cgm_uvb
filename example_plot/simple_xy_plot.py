#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

def plot_xy(x, y, fig_name = 'test.pdf', labels= None, vline = None, hline = None, show = False,
            xrange = None, figure_size = [6, 5]):
    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 11}
    plt.rc('font', **font)
    # matplotlib.rcParams['axes.linewidth'] = 1.5

    fig, (ax1) = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))

    ax1.plot(x, y, marker = '.')

    if vline is not None:
        ax1.axvline(vline, ls = '--')
    if hline is not None:
        ax1.axhline(hline, ls = '--')

    if labels is not None:
        ax1.set_xlabel(labels[0])
        ax1.set_ylabel(labels[1])

    if xrange is not None:
        ax1.set_xlim(xrange[0], xrange[1])

    # decorating the plot
    ax1.tick_params(direction='in', length=5, width=1.5)
    ax1.tick_params(direction='in', which='minor', length=3.5, width=1.5)
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)

    fig.tight_layout(h_pad=1)
    fig.savefig(fig_name, bbox_inches='tight')
    if show :
        plt.show()

"""
# example to use; type following in the ipython terminal
from CGM_UVB.example_plot.simple_xy_plot import plot_xy
import numpy as np
x = np.arange (20)
y = x**2
plot_xy(x=x, y=y, show = True)
"""