# to test mcmc routine for fitting a straight line from points
"""
Resources
http://jakevdp.github.io/blog/2014/06/14/frequentism-and-bayesianism-4-bayesian-in-python/
https://emcee.readthedocs.io/en/stable/tutorials/line/
"""
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner

def generate_data_and_plot( b_true =2, m_true = -1, number_of_points = 20, plot = False, seed =123, same_sigma = None):
    # generating a data on straight line
    np.random.seed(seed)
    x = 10 * np.random.rand (number_of_points) # x array
    y =  m_true * x + b_true  # y from equation
    # error sigmas say 0.5 times a uniform random values between (0.5, 2)
    if same_sigma is None:
        sigma_y = np.repeat(0.5, number_of_points)* np.random.uniform (0.5, 2, number_of_points)
    else:
        sigma_y = same_sigma
    y_noise = sigma_y * np.random.randn(number_of_points) # sigma_y times a gaussian of mean 0 and sigma 1
    y = y + y_noise

    if plot:
        plt.errorbar(x, y, yerr=sigma_y, fmt=".k", capsize=0)
        x0 = np.linspace(0, 10, 500)
        plt.plot(x0, m_true * x0 + b_true, "k", alpha=0.3, lw=3)
        plt.xlim(0, 10)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.show()

    return x, y, sigma_y

#----- for fitting the data we need log posterior i.e product of log likelihood and log prior
def log_likelihood(theta, x, y, yerr):
    """
    For a gaussian distributed errors
    :param theta: parameters
    :param x: data x
    :param y: data y
    :param yerr: data err
    :return:
    """
    b, m = theta # parameters
    y_model = b + m * x # get y value for parameter from model
    lnL = -0.5 * np.sum(np.log(2 * np.pi * yerr ** 2) + (y - y_model) ** 2 / yerr ** 2)

    return lnL

def log_prior(theta):
    b, m =  theta
    # flat prior
    if -5.0 < m < 5 and 0.0 < b < 10.0 :
        return 0.0
    return -np.inf

def log_posterior(theta, x, y, yerr):
    return log_prior(theta) + log_likelihood(theta, x, y, yerr)


#------------------ here is a way to run code
truths = [2, -1] # (b, m) true values
x, y, yerr = generate_data_and_plot(b_true= truths[0], m_true=truths[1], same_sigma= 0.5)
# Here we'll set up the computation. emcee combines multiple "walkers",
# each of which is its own MCMC chain. The number of trace results will
# be nwalkers * nsteps

ndim = 2  # number of parameters in the model
nwalkers = 50  # number of MCMC walkers
nsteps = 4000  # number of MCMC steps to take

# set theta near the maximum likelihood, with
np.random.seed(0)
starting_guesses = 1e-4 * np.random.randn(nwalkers, ndim) # initialise at a tiny sphere

# Here's the function call where all the work happens:
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(x, y, yerr))

sampler.run_mcmc(starting_guesses, nsteps, progress=True)


