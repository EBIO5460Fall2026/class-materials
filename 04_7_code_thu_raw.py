import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


# Poisson distribution emerging from searching for prey
# Number of prey items found in a set time

def num_prey_found(s, t_max, rng):

    t = 0
    prey = 0
    while t < t_max:

        # draw waiting time until next prey found
        w = rng.exponential(scale=1/s)

        # advance time
        t = t + w

        # event occurs
        prey = prey + 1

    return prey - 1


# Test function

num_prey_found(s, t_max, rng)


# Simulation study

rng = np.random.default_rng(20598)

sims = 100000

s = 0.05          # search rate (ca probability per unit time)
t_max = 500       # stop simulation at this time
print(1/s)        # average time to find next prey item

nprey = np.empty(sims)

for i in range(sims):
    nprey[i] = num_prey_found(s, t_max, rng)

# Distribution (PMF) of number found
plt.figure()
bins = np.arange(np.min(nprey)-0.5, np.max(nprey)+0.5+1)
plt.hist(nprey, bins, density=True)
plt.xlabel("Number prey found")
plt.ylabel("Probability mass")

# Add theoretical PMF to plot
X = np.arange(np.min(nprey), np.max(nprey)+1)
pmf = stats.poisson.pmf(X, mu = s * t_max)
plt.plot(X, pmf, marker='o', linestyle='none', fillstyle='none')
plt.title("Poisson distribution (points) vs simulation (bars)")



# Simulating multiple event types
# Prey search and death

rng = np.random.default_rng(1234)

s = 365 / 3   # Catches prey every few days
m = 1 / 15    # Expected lifetime 15 years

alive = True
t = 0
prey = 0
while alive:

    # total rate (intensity)
    intensity = s + m

    # waiting time until next event
    w = rng.exponential(scale=1/intensity)

    # advance time
    t = t + w

    # determine which event occurred
    p_prey = s / intensity

    if rng.binomial(1, p_prey):

        # prey found
        prey = prey + 1

    else:

        # death
        alive = False

print("Time of death:", t)
print("Prey found:", prey)



# Geometric distribution emerging from searching for prey + death
# Number of prey items over a lifetime

def find_prey_then_die(s, m, rng):

    # overall event rate
    intensity = s + m

    # Initial state
    alive = True
    t = 0
    prey = 0

    # Dynamics
    while alive:

        # waiting time until next event
        w = rng.exponential(scale=1/intensity)

        t = t + w

        # determine which event occurred
        p_prey = s / intensity

        if rng.binomial(1, p_prey):

            # prey found
            prey = prey + 1

        else:

            # death
            alive = False

    return prey, t



# Test function

find_prey_then_die(s, m, rng)



# Simulation study

rng = np.random.default_rng(20598)

sims = 100000

# Imagine a spider catching insects
s = 1           # Search (catch) rate (per day)
m = 1 / 14      # Death rate (days)
print(1 / s)    # Expected time to find next prey item
print(1 / m)    # Expected lifetime

nprey = np.empty(sims)
age_at_death = np.empty(sims)

for i in range(sims):
    nprey[i], age_at_death[i] = find_prey_then_die(s, m, rng)
    if i % 100 == 0:
        print(i)

# Distribution (PMF) of number prey found
plt.figure()
bins = np.arange(np.min(nprey)-0.5, np.max(nprey)+0.5+1)
plt.hist(nprey, bins, density=True)
plt.xlabel("Number prey found")
plt.ylabel("Probability mass")

# Add theoretical PMF to plot
X = np.arange(np.min(nprey), np.max(nprey)+1)
# Scipy parameterization of is num trials to first success
# we have number of failures *before* success, so add 1
pmf = stats.geom.pmf(X + 1, m / (s + m))
plt.plot(X, pmf, marker='o', linestyle='none', fillstyle='none')
plt.title("Geometric distribution (points) vs simulation (bars)")

# Expected number of prey E[N] = s / m
print(np.mean(nprey))
print(s / m)

# Variance in number of prey E[N] = s / m
print(np.var(nprey, ddof=1))
print(s * (s + m) / m**2)

