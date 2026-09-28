# Diffusion approximation to the continuous-time birth-death process
# Deterministic skeleton + stochastic deviation
# Now N is a positive continuous random variable

import numpy as np
import matplotlib.pyplot as plt

b = 0.25
m = 0.2
N_0 = 10
dt = 0.01
t_max = 20

rng = np.random.default_rng()
t = 0
N_t = N_0
N = [N_t]
times = [t]
while t <= t_max:
    t = t + dt
    Z = rng.normal() #standard normal deviate
#   N(t + dt) = (deterministic skeleton) + (stochastic deviation)
    N_t = ( N_t + (b - m) * N_t * dt )  +  ( np.sqrt((b + m) * N_t * dt) * Z )
    times.append(t)
    N.append(N_t)

plt.figure()
plt.plot(times, N)
plt.ylim(bottom=0)
plt.xlabel("Time")
plt.ylabel("N")
plt.title("Birth-death stochastic approximation")


# Package into a function

# Including rewriting for some modest (30%) speed gains (known array size,
# pre-allocate a python list, pregenerate a Z array, break if extinct, don't
# track times). Further speed up would be possible using Numba for just-in-time
# compiling. I checked this for exact match to above prototype with same seed.

def birth_death_diff_approx(N_0, b, m, dt, t_max, rng):
    n_steps = int(t_max / dt)
    N = [0.0] * (n_steps + 1)
    N[0] = N_0
    Z = rng.normal(size=n_steps)
    for i in range(n_steps):
        N[i+1] = N[i] + (b - m) * N[i] * dt + np.sqrt((b + m) * N[i] * dt) * Z[i]
        if N[i+1] < 0:
            N[i+1] = 0 #extinct
            break
    return N


# Simulation study
# Does the stochastic approximation give similar mean variance and extinction
# compared to the exact Gillespie IBM?

sims = 1000
b = 0.25
m = 0.2
N_0 = 20
dt = 0.01
t_max = 100

rng = np.random.default_rng(69948)
N = []
for i in range(sims):
    N_i = birth_death_diff_approx(N_0, b, m, dt, t_max, rng)
    N.append(N_i)
    if i % 100 == 0:
        print(i)

n_steps = int(t_max / dt)
t_grid = np.linspace(0, t_max, n_steps + 1)
N = np.array(N)

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
for r in range(25):
    plt.plot(t_grid, N[r, :], lw=0.75, alpha=0.5, color='C0')
plt.plot(t_grid, N_mean, color='red', label='mean')
plt.legend()
plt.ylim(bottom=0)
#plt.ylim(top=7000)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death diffusion approx: N_0={N_0} b={b} m={m}")

# Code test: compare to deterministic model of exponential growth
# Right on
N_det = N_0 * np.exp((b - m) * t_grid) #b - m = 0.05
plt.plot(t_grid, N_det, color='red', linestyle='--', label='deterministic model')
plt.legend()

# Plot mean and standard deviation (process variability)
plt.figure()
plt.fill_between(t_grid, N_mean - sd, N_mean + sd, alpha=0.3, color='C0',
    label='standard deviation')
for r in range(25):
    plt.plot(t_grid, N[r, :], lw=0.75, alpha=0.5)
plt.plot(t_grid, N_mean, color='C0', label='mean')
plt.legend()
plt.ylim(bottom=0)
plt.ylim(top=7000)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death diffusion approx: N_0={N_0} b={b} m={m}")

# Plot extinct populations
plt.figure()
for r in np.where(np.min(N, axis=1) == 0 )[0]:
    plt.plot(t_grid, N[r, :], lw=0.75, alpha=0.5)
plt.plot(t_grid, N_mean, color='C0', label='mean')
plt.legend()
plt.ylim(bottom=0)
plt.ylim(top=60)
plt.xlim(right=t_max)
plt.xlabel("Time")
plt.ylabel("N")
plt.title(f"Birth-death diffusion approx: N_0={N_0} b={b} m={m}")
