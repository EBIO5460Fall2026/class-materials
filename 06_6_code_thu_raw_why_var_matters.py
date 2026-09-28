# Continuous-time birth-death IBM
# Simulation study: mean and variance of the ensemble over time and extinction.
# Shows that the mean model misses important biological phenomena.
# Scenarios with the same mean outcome can have very different variance
# and extinction risk.

import numpy as np
import matplotlib.pyplot as plt

# Function to grow the size of arrays
# Factor=2 (doubling size) by default
# Suggested fills for various dtypes
#     dtype    fill
#    ----------------------
#     float    0 or np.nan
#     int      -1
#     bool     False

def grow(x, fill, factor=2):

    # get original array size
    n = len(x)

    # make new bigger array (double by default)
    x_new = np.full(factor * n, fill, dtype=x.dtype)

    # copy original array into new one
    x_new[:n] = x

    return x_new


def birth_death_ibm(alive, b, m, t_max, rng):

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
    m_intensity = np.sum(m * alive)
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
            m_intensity = m_intensity + m[i]
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
            m_intensity = m_intensity - m[i]
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


# Simulation study
# Mean and variance of the ensemble over time and extinction
# Investigation: change parameters and initial state

# Simulation control
rng = np.random.default_rng(29732)
sims = 1000
t_max = 100

# Parameters and initial state
N_0 = 25
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate

# Storage for result
result = []

# Realizations
for i in range(sims):
    N_series, t = birth_death_ibm(alive, b, m, t_max, rng)
    result.append((N_series, t))
    if i % 100 == 0:
        print(i)

# Compile results onto a grid of times
# rows = realizations, columns = grid times
t_grid = np.linspace(0, 100)
N = np.empty((sims, len(t_grid)), dtype=int)

for r in range(sims):
    N_series, t = result[r]
    for c in range(len(t_grid)):

        # Find first event time after the grid time
        # unless it's the last time in the series
        j = 0
        while t[j] <= t_grid[c] and j < (len(t) - 1):
            j += 1

        # Record N from the final or previous event
        if j != (len(t) - 1):
            j = j - 1
        N[r, c] = N_series[j]

    if r % 100 == 0:
        print(r)


# Mean, standard deviation of realizations, standard error of mean
N_mean = np.mean(N, axis=0)
sd = np.std(N, axis=0, ddof=1)
se = sd / np.sqrt(sims)

# Number extinct
extinct = sum(np.min(N, axis=1) == 0)
print(f"Number extinct: {extinct}")
print(f"Proportion extinct: {extinct / sims}")

# Plot mean and standard error
plt.figure()
plt.fill_between(t_grid, N_mean - 2 * se, N_mean + 2 * se, alpha=0.3,
    color='red', label='s.e.')
# first 25 realizations (from original results list, not regular grid)
for i in range(25):
    N_step, t = result[i]
    plt.step(t, N_step, where='post', lw=0.75, alpha=0.5, color='C0')
plt.plot(t_grid, N_mean, color='red', label='mean')
plt.legend()
plt.ylim(bottom=0)
#plt.ylim(top=7000)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death IBM: N_0={N_0} b={b[0]} m={m[0]}")

# Compare to deterministic model of exponential growth
# The deterministic model emerges as the mean of the stochastic process
N_det = N_0 * np.exp((b[0] - m[0]) * t_grid) #b - m = 0.05
plt.plot(t_grid, N_det, color='black', linestyle='--', label='deterministic')
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
plt.ylim(top=7000)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death IBM: N_0={N_0} b={b[0]} m={m[0]}")

# Plot extinct populations
plt.figure()
for i in np.where(np.min(N, axis=1) == 0 )[0]:
    N_step, t = result[i]
    plt.step(t, N_step, where='post', lw=0.75, alpha=0.5)
plt.plot(t_grid, N_mean, color='C0', label='mean')
plt.legend()
plt.ylim(bottom=0)
plt.ylim(top=60)
plt.xlim(right=t_max)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death IBM: N_0={N_0} b={b[0]} m={m[0]}")
