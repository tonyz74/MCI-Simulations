import scipy.integrate as integrate
import numpy as np
from scipy.stats import norm

# Define the bounds
lower_bound_x0, upper_bound_x0 = -3, 4
lower_bound_x1, upper_bound_x1 = -1, 2

# Define the Gaussian PDF for the region
def gaussian_pdf(x0, x1, sigma):
    return norm.pdf(x0, scale=sigma) * norm.pdf(x1, scale=sigma)

# Integrate the PDF over the given bounds
sigma = 0.5
NC, _ = integrate.dblquad(lambda x1, x0: gaussian_pdf(x0, x1, sigma),
                          lower_bound_x0, upper_bound_x0,
                          lower_bound_x1, upper_bound_x1)

print("Normalization constant (NC):", NC)
