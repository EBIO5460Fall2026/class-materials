# Approximate Bayesian Computation
# Demonstrate the basic rejection ABC algorithm
# Train the logistic IBM

import numpy as np
import matplotlib.pyplot as plt

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


# Simulate a single data set
# This will be the data we wish to analyze.
# Our goal will be to construct the posterior distribution for the parameter c.
# For simplicity we'll assume N_0, b, m are known, while only c is unknown, but
# the general algorithm here extends easily to multiple unknown parameters.

# Simulation control
rng = np.random.default_rng(29732)
t_max = 200

# Parameters and initial state
N_0 = 25
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate
c_true = 0.001              # per capita competition effect on mortality

# Simulate IBM and snap to equal spaced time grid
N_obs, t_obs = birth_death_ibm_dd(alive, b, m, c_true, t_max, rng)

# Observed summary statistic (e.g. the data itself on a regular grid)
t_grid = np.linspace(0, 200)   # length 50 vector
N_obs_grid = regularize_series(N_obs, t_obs, t_grid)

# Here's what it looks like
plt.figure()
plt.plot(t_grid, N_obs_grid)
plt.xlabel("Time")
plt.ylabel("N")
plt.ylim(bottom=0)


# ABC rejection algorithm
# Basic algorithm is 7 lines of code
# Train c

n_samples = 10000
tol = 15 #roughly top 1% of simulations
c_posterior = []
for i in range(n_samples):

    # draw c from prior
    c = rng.uniform(0.0002, 0.003)

    # simulate the model
    N_sim, t_sim = birth_death_ibm_dd(alive, b, m, c, t_max, rng)

    # summary statistic
    N_sim_grid = regularize_series(N_sim, t_sim, t_grid)

    # distance between statistics
    d = np.sqrt(np.mean((N_sim_grid - N_obs_grid)**2))

    # accept or reject as posterior sample
    if d < tol:
        c_posterior.append(c)

    # monitor progress
    if i % 1000 == 0:
        print(i)


# Posterior (and prior) distribution
c_prior = rng.uniform(0.0002, 0.003, n_samples)
plt.figure()
plt.hist(c_prior, bins=40, density=True, alpha=0.4, label="prior")
plt.hist(c_posterior, bins=20, density=True, alpha=0.7, label="posterior")
plt.axvline(c_true, color="blue", label="true value")
plt.xlabel("c")
plt.ylabel("Density")
plt.legend()

print(len(c_posterior))      # number accepted
print(np.mean(c_posterior))  # posterior mean
print(c_true)                # true c


# Here is a Visualization of the rejection algorithm
# 500 iterations

from matplotlib.lines import Line2D

plt.figure()
for i in range(500):

    # draw c from prior
    c = rng.uniform(0.0002, 0.003)

    # simulate the model
    N_sim, t_sim = birth_death_ibm_dd(alive, b, m, c, t_max, rng)

    # summary statistic
    N_sim_grid = regularize_series(N_sim, t_sim, t_grid)

    # distance between statistics
    d = np.sqrt(np.mean((N_sim_grid - N_obs_grid)**2))

    # accept or reject as posterior sample
    if d < tol:
        plt.plot(t_grid, N_sim_grid, color="C1", alpha=0.5)
    else:
        plt.plot(t_grid, N_sim_grid, color="lightgray", alpha=0.1)

plt.plot(t_grid, N_obs_grid, color="black")

legend_elements = [
    Line2D([0], [0], color="black", label="Observed"),
    Line2D([0], [0], color="C1", label="Accepted"),
    Line2D([0], [0], color="lightgray", label="Rejected")
]

plt.legend(handles=legend_elements)
plt.xlabel("Time")
plt.ylabel("N")
