
import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import emcee
import corner


#----data
def get_true_model(model_Q, Q= 17):
    """
    :param model: The data where Q18 model is stored
    :return: a row of ion column densities at n_H = 1e-4 cm^-2
    """
    #model = model_Q.split('_Q')[0] + '_Q{}.fits'.format(Q)
    model = '/home/vikram/cloudy_run/hst/NH15/Q17/Q{}.fits'.format(Q)
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


#----for mcmc
def log_likelihood(theta, interp_logf, obs_ion_col, col_err, reference_log_metal = -1.0):
    """
    For a gaussian distributed errors
    :param theta: parameters [nH, Z]
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    lognH, logZ =  theta
    # get metal ion column density for n_H and Z = 0.1
    col = 10 ** interp_logf(lognH)
    # scale the column densities by the metallicity Z
    metal_scaling_linear = 10 ** logZ / 10 ** reference_log_metal
    model_col = np.log10(col * metal_scaling_linear)

    lnL = -0.5 * np.sum(np.log(2 * np.pi * col_err ** 2) + (obs_ion_col - model_col) ** 2 / col_err ** 2)

    return lnL

def log_prior(theta):
    lognH, logZ =  theta
    # flat prior
    if -6 < lognH < -2 and -2 < logZ < 1 :
        return 0.0
    return -np.inf

def log_posterior(theta, interp_func, data_col, sigma_col):
    log_p = log_prior(theta) + \
            log_likelihood(theta, interp_logf = interp_func, obs_ion_col = data_col, col_err = sigma_col)

    return log_p


def run_mcmc(model_Q, ions_to_use, true_Q =18, figname = 'test.pdf', same_error = True, sigma_noise = 0.1):
    # run_mcmc(model_Q= model, ions_to_use= ions)
    # ------------------ here is a way to run code
    truths = [-4, -1]  # (lognH, logZ) true values
    number_of_ions = len(ions_to_use)


    #np.random.seed(0)
    if same_error:
        sigma_col = sigma_noise * np.ones(number_of_ions)
        noise  = np.random.normal(0, sigma_noise, number_of_ions)
        print(noise)
    else:
        sigma_col = np.random.uniform(0.1, 0.3, number_of_ions)


    data_col_all = get_true_model(model_Q, Q=true_Q)

    # converting astropy table row to a list
    data_col = []
    for name in ions_to_use:
        data_col.append(data_col_all[name][0])

    print(data_col)
    # shifting the column values
    data_col =  10**(np.log10(data_col) + noise)
    print(data_col)

    print(np.log10(data_col), sigma_col)

    interp_logf = get_interp_func(model_Q, ions_to_use)

    # Here we'll set up the computation. emcee combines multiple "walkers",
    # each of which is its own MCMC chain. The number of trace results will
    # be nwalkers * nsteps

    ndim = 2  # number of parameters in the model
    nwalkers = 50  # number of MCMC walkers
    nsteps = 5000  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    n_guess = np.random.uniform(-5, -2, nwalkers)
    z_guess = np.random.uniform(-2, 1, nwalkers)
    starting_guesses = np.vstack((n_guess, z_guess)).T  # initialise at a tiny sphere

    # Here's the function call where all the work happens:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(interp_logf, np.log10(data_col), sigma_col))
    sampler.run_mcmc(starting_guesses, nsteps, progress=True)

    # find out number of steps
    tau = sampler.get_autocorr_time()  # number of steps needed to forget the starting position
    #print(tau)
    thin = int(np.mean(tau) / 2)  # use this number for flattning the sample as done below
    #thin = 100
    flat_samples = sampler.get_chain(discard=thin * 20, thin= 5, flat=True)
    # we are discarding some initial steps roughly 5 times the autocorr_time steps
    # then we thin by about half the autocorrelation time steps for plotting => one does not have to do this step

    labels = ['log nH', 'log Z']
    uvb_q= 17 #int((model_Q.split('try_Q')[-1]).split('.fits')[0])

    if uvb_q == true_Q:
        fig = corner.corner(flat_samples, labels=labels, truths=truths, quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12})
    else:
        fig = corner.corner(flat_samples, labels=labels, quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_kwargs={"fontsize": 12})

    fig.savefig(figname)

    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        print(labels[i], '=', mcmc[1], q[0], q[1])


    return flat_samples, ndim


ions_to_use= ['C+3', 'N+3', 'Si+3', 'O+5', 'C+2']
true_Q =17
outpath = '/home/vikram/cloudy_run/hst'

outfile = outpath + '/NH15.fits'



uvb_array= [17]
out_tab =  tab.Table()
for uvb_q in uvb_array:
    model_Q = '/home/vikram/cloudy_run/hst/NH15/Q17/Q{}.fits'.format(uvb_q)
    name = 'Q17'
    figname = outpath + '/' + name + '.pdf'

    flat_samples, ndim = run_mcmc(model_Q=model_Q, ions_to_use=ions_to_use, true_Q=true_Q, figname=figname)
    # to efficiently save numpy array
    save_file_name = outpath + '/' + name
    np.save(save_file_name, flat_samples)

    out =[[uvb_q]]
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        out.append([mcmc[1]])
        out.append([q[0]])
        out.append([q[1]])

    print(out)
    t = tab.Table(out, names = ('Q', 'nH', 'n16', 'n84', 'Z', 'Z16', 'Z84'))
    out_tab = tab.vstack((out_tab, t))


out_tab.write(outfile, overwrite = True)

