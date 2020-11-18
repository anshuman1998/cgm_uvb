# to test mcmc routine for fitting a straight line from points
import numpy as np
import matplotlib.pyplot as plt

def generate_data_and_plot(m = -1, b =2, number_of_points = 20, plot = False):
    x = 10 * np.random.rand (number_of_points) # x array
    y =  m * x + b  # y from equation
    # error sigmas say 0.2 times a uniform random values between (0.5, 2)
    sigma_y = np.repeat(0.2, number_of_points)* np.random.uniform (0.5, 2, number_of_points)
    y_noise = sigma_y * np.random.randn(number_of_points) # sigma_y times a gaussian of mean 0 and sigma 1
    y = y + y_noise

    if plot:


