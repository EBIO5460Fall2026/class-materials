# Simulation: continuous time algorithm

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom

# Initiate the random number generator

rng = np.random.default_rng(2293)

m = 0.02         # roughly probability of death per year (so, lives ca 50 years)
delta_t = 1/365  # time increment 1 day, units are years
p = m * delta_t  # probability per day
print(p)

state = 0
t = 0
while state == 0:
    t = t + delta_t
    if ( rng.binomial(1, p, 1) ):
        state = 1
print(t)
