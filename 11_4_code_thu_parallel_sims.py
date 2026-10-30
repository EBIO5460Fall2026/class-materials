# Create a large collection of labelled simulations using parallel computation

import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor #parallel library
from tqdm import tqdm #for tracking progress
import os

# How many cores?
print(os.cpu_count())

# Helper function to grow the size of arrays
#
def grow(x, fill, factor=2):
    n = len(x)
    x_new = np.full(factor * n, fill, dtype=x.dtype)
    x_new[:n] = x
    return x_new

# Logistic IBM (from week 7)
# Continuous-time birth-death with density dependent mortality
# m_j = m + cN
# c: increase in death rate caused by each individual
#
def birth_death_ibm_dd(alive, b, m, c, t_max, rng):

    # Initial biological state

    t = 0.0      #needs to be float, e.g 0.0

    # --population scale (automatic)
    N_t = np.sum(alive)     #scale transition (sum over individuals -> population)
    N = np.array([N_t])     #track N
    times = np.array([t])   #track time

    #  --keep a list of which individuals are alive
    alive_list = np.arange(N_t)

    # --calculate initial event intensity
    b_intensity = np.sum(b * alive)
    m_intensity = np.sum(m * alive + c * N_t)                              #<new
    intensity = b_intensity + m_intensity

    # Grow arrays to make room for new data (ca 6% time gain)

    alive = grow(alive, fill=False, factor=100)
    b = grow(b, fill=0, factor=100)
    m = grow(m, fill=0, factor=100)

    N = grow(N, fill=-1, factor=200)
    times = grow(times, fill=np.nan, factor=200)

    alive_list = grow(alive_list, fill=-1, factor=10)

    # Simulate events

    isize = len(alive)
    psize = len(N)
    irow = N_t   #array row we insert the next individual info into
    prow = 1     #array row we insert the next population info into

    while t < t_max and N_t > 0:

        # grow arrays if needed
        if irow == isize:
            alive = grow(alive, 0)
            b = grow(b, 0)
            m = grow(m, 0)
            isize = isize * 2

        if prow == psize:
            N = grow(N, -1)
            times = grow(times, np.nan)
            psize = psize * 2

        if len(alive_list) == N_t:
            alive_list = grow(alive_list, -1)

        # waiting time to next event
        w = rng.exponential(scale=1/intensity)

        # advance time
        t = t + w

        # which event type occurs?
        event_birth = rng.binomial(1, b_intensity / intensity) #1 = birth

        # which individual does the event happen to?
        j = rng.integers(N_t)
        i = alive_list[j]

        # update system state (do event)
        if event_birth:
            # insert new individual in arrays
            alive[irow] = 1
            b[irow] = b[i]  #inherit from parent
            m[irow] = m[i]
            # update alive list
            alive_list[N_t] = irow
            # update intensities
            b_intensity = b_intensity + b[i]
            m_intensity = m_intensity + m[i] + c * (2 * N_t + 1)           #<new
            # update N
            N_t = N_t + 1
            # increment row index for next time
            irow += 1
        else:
            # death
            alive[i] = 0
            # update alive list (swap with last)
            alive_list[j] = alive_list[N_t - 1]
            # update intensities
            b_intensity = b_intensity - b[i]
            m_intensity = m_intensity - m[i] - c * (2 * N_t - 1)           #<new
            # update N
            N_t = N_t - 1

        intensity = b_intensity + m_intensity

        # insert population level states in arrays
        N[prow] = N_t
        times[prow] = t
        # increment row index for next time
        prow += 1

    # Trim arrays (remove the unused elements)
    # alive = alive[:irow]
    # b = b[:irow]
    # m = m[:irow]

    if t > t_max:
        prow = prow - 1

    N = N[:prow]
    times = times[:prow]

    return N, times

# Helper function to put simulation data on a regular (even spaced) time grid
# Processes a single realization from the IBM
#
def regularize_series(N_series, times, t_grid):
    N_grid = np.empty(len(t_grid), dtype=int)
    for col in range(len(t_grid)):
        j = 0
        while times[j] <= t_grid[col] and j < (len(times) - 1):
            j += 1
        if j != (len(times) - 1):
            j -= 1
        N_grid[col] = N_series[j]
    return N_grid


# A convenient approach is to start by creating a large collection of
# simulations using parallel computation

# Globals
t_max = 200
N_0 = 25
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate

# Package one simulation in a function
# Each call to the function will get its own seed
#
def abc_draw(seed):

    rng = np.random.default_rng(seed)
    c = rng.uniform(0.0002, 0.003)
    N_sim, t_sim = birth_death_ibm_dd(alive, b, m, c, t_max, rng)
    N_sim_grid = regularize_series(N_sim, t_sim, t_grid)

    return c, N_sim_grid

# In this approach, each worker gets a copy of the global environment, so it has
# all the functions and variables that have been initialized. The map() function
# does the work iterating over seeds. It's functional programming equivalent to:
#     for s in seeds:
#         results.append(abc_draw(s))
# tqdm() adds a progress tracker.

n_draws = 100000
main_rng = np.random.default_rng(40983) #set this for reproducibility
seeds = main_rng.integers(0, 2**32, size=n_draws)
with ProcessPoolExecutor(max_workers=14) as pool:
    results = list(
        tqdm(
            pool.map(abc_draw, seeds)
        )
    )

# Pythonic way to reorganize results (* unpacks, zip groups)
c_sample, N_sim_grid = zip(*results)
c_sample = np.array(c_sample)
N_sim_grid = np.array(N_sim_grid)

# Save results to file
np.save("c_sample.npy", c_sample)
np.save("N_sim_grid.npy", N_sim_grid)

# While we're at it, save the true dataset to file also
rng = np.random.default_rng(29732)
c_true = 0.001
N_obs, t_obs = birth_death_ibm_dd(alive, b, m, c_true, t_max, rng)
t_grid = np.linspace(0, 200)
N_obs_grid = regularize_series(N_obs, t_obs, t_grid)
np.save("N_obs_grid.npy", N_obs_grid)

