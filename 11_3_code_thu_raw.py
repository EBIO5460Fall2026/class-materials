# Approximate Bayesian Computation
# Demonstrate the basic algorithm
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

        # insert population

# Simulation study
# Mean and variance of the ensemble over time and extinction
# Investigation: change parameters and initial state

# Simulation control
rng = np.random.default_rng(29732)
sims = 1000
t_max = 200

# Parameters and initial state
N_0 = 25
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate
c = 0.001                  # per capita competition effect on mortality

# Storage for result
result = []

# Realizations
for i in range(sims):
    N_series, t = birth_death_ibm_dd(alive, b, m, c, t_max, rng)
    result.append((N_series, t))
    if i % 100 == 0:
        print(i)

# Compile results onto a grid of times
# rows = realizations, columns = grid times
t_grid = np.linspace(0, 200)


# Mean, standard deviation of realizations, standard error of mean
N_mean = np.mean(N, axis=0)
sd = np.std(N, axis=0, ddof=1)
se = sd / np.sqrt(sims)

# Mean, standard deviation of realizations, standard error of mean
# conditioned on being extant
N_extant = N[np.min(N, axis=1) != 0]
N_mean_extant = np.mean(N_extant, axis=0)
sd_extant = np.std(N_extant, axis=0, ddof=1)
se_extant = sd_extant / np.sqrt(len(N_extant))

# Number extinct
extinct = sum(np.min(N, axis=1) == 0)
print(f"Number extinct: {extinct}")
print(f"Proportion extinct: {extinct / sims}")


# Plot means and standard errors
plt.figure()
plt.fill_between(t_grid, N_mean - 2 * se, N_mean + 2 * se, alpha=0.3,
    color='red', label='s.e.')
plt.fill_between(t_grid, N_mean_extant - 2 * se_extant,
    N_mean_extant + 2 * se_extant, alpha=0.3, color='C0', label='s.e. extant')
# first 25 realizations (from original results list, not regular grid)
for i in range(10):
    N_step, t = result[i]
    plt.step(t, N_step, where='post', lw=0.75, alpha=0.5)
plt.plot(t_grid, N_mean, color='red', label='mean')
plt.plot(t_grid, N_mean_extant, color='C0', label='mean extant')
plt.legend()
plt.ylim(bottom=0)
#plt.ylim(top=7000)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death IBM: N_0={N_0} b={b[0]} m={m[0]} c={c}")


# Compare to deterministic model of logistic growth
# The deterministic model emerges as the mean of the stochastic process
K = (b[0] - m[0]) / c
N_det = N_0 * K / (N_0 + (K - N_0) * np.exp(-1 * (b[0] - m[0]) * t_grid))
plt.plot(t_grid, N_det, color='black', linestyle='--', label='deterministic')
plt.legend()

# Expected equilibrium
E_N = K - np.mean(sd_extant[25:49]**2) / K
plt.plot(np.arange(101,201), np.full(100, E_N), color='purple',
    linestyle='--', label='expected')
plt.legend()


# Plot mean and standard deviation (process variability)
plt.figure()
plt.fill_between(t_grid, N_mean - sd, N_mean + sd, alpha=0.3, color='C0',
    label='standard deviation')
for i in range(25):
    #plt.plot(t_grid, N[i, :], lw=0.75, alpha=0.5) #alt: grid times
    N_step, t = result[i]
    plt.step(t, N_step, where='post', lw=0.75, alpha=0.3)
plt.plot(t_grid, N_mean, color='C0', label='mean')
plt.legend()
plt.ylim(bottom=0)
plt.ylim(top=100)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death IBM: N_0={N_0} b={b[0]} m={m[0]} c={c}")


# Plot extinct populations
plt.figure()
for i in np.where(np.min(N, axis=1) == 0 )[0]:
    N_step, t = result[i]
    plt.step(t, N_step, where='post', lw=0.75, alpha=0.5)
plt.plot(t_grid, N_mean, color='C0', label='mean')
plt.legend()
plt.ylim(bottom=0)
plt.ylim(top=100)
plt.xlim(right=t_max)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death IBM: N_0={N_0} b={b[0]} m={m[0]} c={c}")
 level states in arrays
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


# Helper function to put a single realization on a regular grid
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


# ------------------------------------------------------------------------------
# Simulated data set

rng = np.random.default_rng(33093)
t_max = 200

# "True" parameters
N_0 = 25
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate
c_true = 0.001              # per capita competition effect on mortality

# Generate data
N_obs, t_obs = birth_death_ibm_dd(alive, b, m, c_true, t_max, rng)

# Observed summary statistic (e.g. the data itself on a regular grid)
t_grid = np.linspace(0, 200, 50)
N_obs_grid = regularize_series(N_obs, t_obs, t_grid)
