
import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.table as tab
import matplotlib.pyplot as plt


#----data
def get_true_model(model_Q, Q= 18):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    model = model_Q.split('_Q')[0] + '_Q{}.fits'.format(Q)
    data = tab.Table.read(model)
    true_ion_col = data [data['hden'] == 1e-4]
   # print(true_ion_col)

    return true_ion_col

#----model interpolation
def get_interp_func(model_Q, ions_to_use):
    number_of_ions = len(ions_to_use)

    model = tab.Table.read(model_Q)
    sorted_model = model[ions_to_use]
    hden_array = np.array(model['hden'])

    model_touple = ()
    for j in range(number_of_ions):
        model_touple += (sorted_model[ions_to_use[j]],)

    # interpolating in log log scale
    logf = interp1d(np.log10(hden_array), np.log10(model_touple), fill_value='extrapolate')

    return logf

def sample_plot(npy_file, interp_logf, true_ions, true_col, true_sigma, figname ='test.png', samples = 100,
                reference_log_metal = -1):

    # setting the figure
    font = {'family': 'serif', 'weight': 'normal', 'size': 14}
    plt.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 1.5

    figure_size = [7, 6]
    fig, ax = plt.subplots(1, 1, figsize=(figure_size[0], figure_size[1]))

    flat_samples = np.load(npy_file)
    inds = np.random.randint(len(flat_samples), size=samples)
    for ind in inds:
        sample = flat_samples[ind]
        lognH, logZ = sample
        col = 10 ** interp_logf(lognH)
        # scale the column densities by the metallicity Z
        metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
        model_col = col * metal_scaling_linear
        ax.plot(true_ions, model_col,  "C1", alpha=0.2)

    ax.errorbar(true_ions, true_col, yerr = (true_col*true_sigma)*2.303, color = 'k', elinewidth =2, marker = '.',  markersize = 15,
                linestyle = '', zorder =1000)

    ax.set_xlabel('Ion')
    ax.set_xlabel(r'N$_{\rm ion}$')
    # ax.set_xscale('log')
    ax.set_yscale('log')

    # deco
    ax.tick_params(direction='in', length=5, width=1.5)
    ax.tick_params(direction='in', which='minor', length=3.5, width=1.5)
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    # decorating the plot
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)

    fig.savefig(figname, bbox_inches='tight')




def get_mcmc_sample(model_Q, ions_to_use, npyfile,  true_Q =18, figname = 'test.png', same_error = False):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    number_of_ions = len(ions_to_use)

    data_col_all = get_true_model(model_Q, Q=true_Q)
    # converting astropy table row to a list
    data_col = []
    for name in ions_to_use:
        data_col.append(data_col_all[name][0])


    np.random.seed(0)
    if same_error:
        sigma_col = 0.2 * np.ones(number_of_ions)
    else:
        sigma_col = np.random.uniform(0.1, 0.3, number_of_ions)

    print(np.log10(data_col), sigma_col)

    interp_logf = get_interp_func(model_Q, list(ions_to_use))

    sample_plot(npy_file= npyfile, interp_logf= interp_logf, true_ions= ions_to_use, true_col= data_col,
                true_sigma= sigma_col, figname = figname)

    return

uvb_q= 18
model_Q = '/home/vikram/cloudy_run/anshuman/try_Q{}.fits'.format(uvb_q)
ions_to_use= np.array(['C+3', 'N+3', 'Si+3', 'O+5', 'C+2'])
true_Q =18
#-----------imp sort ions
data_col_all = get_true_model(model_Q, Q=true_Q)
# converting astropy table row to a list
data_col = []
for name in ions_to_use:
    data_col.append(data_col_all[name][0])

ions_to_use = ions_to_use[np.argsort(data_col)]
print(ions_to_use, ': sorted')

npyfile =  '/home/vikram/cloudy_run/figures/anshuman_try_Q{}.npy'.format(uvb_q)

get_mcmc_sample(model_Q, ions_to_use, npyfile, figname= 'test18_sample.png')
