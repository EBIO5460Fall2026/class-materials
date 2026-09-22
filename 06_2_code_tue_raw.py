# Continuous-time birth-death IBM

import numpy as np
import matplotlib.pyplot as plt
from time import time

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

    # Grow arrays to make room for new data

    alive = grow(alive, fill=False, factor=100)
    b = grow(b, fill=0, factor=100)
    m = grow(m, fill=0, factor=100)
    isize = len(alive)

    N = grow(N, fill=-1, factor=200)
    times = grow(times, fill=np.nan, factor=200)
    psize = len(N)

    alive_list = grow(alive_list, fill=-1, factor=10)

    # Simulate events

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
    alive = alive[:irow]
    b = b[:irow]
    m = m[:irow]

    if t > t_max:
        prow = prow - 1

    N = N[:prow]
    times = times[:prow]

    return N, times

# --individual scale (set state and parameters)
N_0 = 10
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate
t_max = 30
rng = np.random.default_rng(95756)
N, t = birth_death_ibm(alive, b, m, t_max, rng)

# Check against original
plt.figure()
plt.step(t, N, where='post')
plt.ylim(bottom=0)
plt.xlabel("Time")
plt.ylabel("N")
plt.title("Birth-death IBM")


# Timing study
# e.g. distribution of N at time 30

rng = np.random.default_rng(95756)
sims = 10
t_max = 30
N_0 = 4000
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate

start_time = time() #timer
for i in range(sims):
    N, t = birth_death_ibm(alive, b, m, t_max, rng)
    #print(i)
end_time = time() #stop timer
print(end_time - start_time)



# Simulation study 1
# e.g. distribution of N at time 30

# Simulation control
rng = np.random.default_rng(84923)

sims = 10000
t_max = 30

# Parameters and initial state
N_0 = 10
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate

# Storage for result
N = np.full(sims, -1)

for i in range(sims):
    N_series, t = birth_death_ibm(alive, b, m, t_max, rng)
    N[i] = N_series[len(N_series)-1]
    if i % 1000 == 0:
        print(i)


# Distribution (PMF) of N at t = 30
plt.figure()
bins = np.arange(np.min(N)-0.5, np.max(N)+0.5+1)
plt.hist(N, bins, density=True)
plt.xlabel("N")
plt.ylabel("Probability mass")

# Mean N at t = 30
mean_N = np.mean(N)
print(f"Observed: {mean_N}")

# MC error, s.e.
mcse = np.std(N, ddof=1) / np.sqrt(sims)
print(f"MCSE: {mcse}")

# Approx 95% CI
print(f"95% CI {mean_N - 2 * mcse:.2f}, {mean_N + 2 * mcse:.2f}")

# Expected N is deterministic exponential growth:
# E[N] = N(0)e^{(b-m)t}
print(f"Expected: {N_0 * np.exp(0.05 * t_max)}")

# Probability of extinction, about 6%
p_extinct = np.sum(N == 0) / sims
print(f"Pr extinction: {p_extinct}")

# MC error, s.e.
mcse = np.std(N == 0, ddof=1) / np.sqrt(sims)
print(f"MCSE: {mcse}")

# Approx 95% CI
print(f"95% CI {p_extinct - 2 * mcse:.4f}, {p_extinct + 2 * mcse:.4f}")




# Simulation study 2
# Mean and variance of the ensemble over time

# Simulation control
rng = np.random.default_rng(9956)

sims = 1000
t_max = 100

# Parameters and initial state
N_0 = 25
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate

# Storage for result
result = []

for i in range(sims):
    N_series, t = birth_death_ibm(alive, b, m, t_max, rng)
    result.append((N_series, t))
    if i % 100 == 0:
        print(i)



# Compile results onto a grid of times.
# Algorithm:
# for each realization
#     for each time on the grid
#         find the first event time after the grid time
#         move back one event (unless it's the final event)
#         record N


# Rows = realizations
# Columns = grid times
N = np.empty((sims, len(t_grid)), dtype=int)

for r in range(sims):

    N_series, t = result[r]

    for k in range(len(t_grid)):

        # Find first event time after the grid time
        # unless it's the last time in the series
        j = 0
        while t[j] <= t_grid[k] and j < (len(t) - 1):
            j += 1

        # Record N from the final or previous event
        if j != (len(t) - 1):
            j = j - 1
        N[r, k] = N_series[j]


N_mean = np.mean(N, axis=0)
se = np.std(N, axis=0, ddof=1) / np.sqrt(sims)


# Now we can plot realizations and the mean with 95% CI

plt.figure()

for r in range(25):
    plt.plot(t_grid, N[r, :], lw=0.75, alpha=0.5)

plt.ylim(bottom=0)
plt.xlabel("Time")
plt.ylabel("N")
plt.title("Birth-death IBM")

plt.plot(t_grid, N_mean, color='C0')

plt.fill_between(t_grid, N_mean - 2 * se, N_mean + 2 * se, alpha=0.3, color='C0')

# Compare to deterministic model of exponential growth
# We see that the det model emerges as the mean of the stochastic process
N_det = N_0 * np.exp(0.05 * t_grid)
plt.plot(t_grid, N_mean, color='red', linestyle='--')
