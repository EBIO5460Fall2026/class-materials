import numpy as np
import matplotlib.pyplot as plt

# RNG from last time

def rand_unif(n, seed):
    a = 16807
    m = 2147483647
    U = np.empty(n)

    I = seed
    for j in range(n):
        I = (a * I) % m
        U[j] = I / m

    return U, I          #also return generator state I


# Start the RNG (same seed as last time)

U, I_state = rand_unif(1, seed=74732)


# Here was our agent-based model from last time

# 100 individuals
# Probability of death in a year = 0.2
# How many die in the first year?

N = 100
alive = np.full(N, 1)    #alive = 1
p_death = np.full(N, 0.2)

# Simulate the death process for each individual

for i in range(N):
    U, I_state = rand_unif(1, I_state)
    if U <= p_death[i]:
        alive[i] = 0    #individual dies

print(sum(alive))



# To investigate this stochastic process, we simulate the agent-based model
# multiple times.

# First, we'll package the individual-based model into a function.

def death_IBM(N, p_d, I_state):

    alive = np.full(N, 1)        #alive = 1
    p_death = np.full(N, p_d)

    for i in range(N):
        U, I_state = rand_unif(1, I_state)
        if U <= p_death[i]:
            alive[i] = 0     #individual dies

    return alive, I_state

# Check that it works (aka one realization of the stochastic process)

alive, I_state = death_IBM(N=100, p_d=0.2, I_state=I_state)
print(alive)
print(sum(alive))


# Now simulate the IBM multiple times

U, I_state = rand_unif(1, seed=924505)

reps = 100000
N = 100
p_d = 0.2
alive = np.empty((N, reps), dtype=int)

for i in range(reps):
    alive[:,i], I_state = death_IBM(N, p_d, I_state)
    # if i % 1000 == 0:  # uncomment for monitoring
    #     print(i)

print(alive)   #each simulation is one column

# Each individual has its own Bernoulli distribution
# Extract all the outcomes for individual 37 (which is row index 37)
# Individuals are indexed 0 - 99 (Python; offset indexing)

individual = 37
dead_int = (alive[individual, :] -1) * -1 #convert alive to dead

plt.figure()
plt.hist(dead_int, bins=[-0.5, 0.5, 1.5], density=True)
plt.xticks([0, 1])
plt.xlabel("Dead")
plt.ylabel("Relative frequency")
plt.show()

# Compare to theoretical PMF
plt.plot([0, 1], [1-p_d, p_d], marker='o', linestyle='none', fillstyle='none')
plt.title("Bernoulli distribution (points) vs simulation (bars)")
plt.show()


# Sum the columns to get total alive or dead in the population

N_alive = np.sum(alive, axis=0)
N_dead = np.sum(np.logical_not(alive), axis=0)
print(N_alive)
print(N_dead)

# The total alive or dead at the population scale is in a sense an emergent
# outcome of the individual-based process. It has its own distribution, which is
# a PMF since it is a discrete (integer) random variable.

# Distribution (PMF) of number dead

plt.figure()
bins = np.arange(np.min(N_dead)-0.5, np.max(N_dead)+0.5+1)
plt.hist(N_dead, bins, density=True)
plt.xlabel("Number dead")
plt.ylabel("Probability mass")
plt.show()

# The distribution of an outcome random variable from an individual-based
# simulation like this can have any arbitrary distribution that may or may not
# be among the standard known or even weird but known distributions. So, the
# individual-based model is a completely general way of obtaining the
# distribution (represented by the histogram). But it's worth being on the
# lookout for known outcomes as it can potentially reduce the complexity of your
# model to use the emergent distribution directly.


# Binomial distribution

# In this case, the outcome for N_dead (or N_alive) is a binomial distribution.
# The binomial distribution arises as the sum of multiple independent Bernoulli
# trials with the same probability. That's exactly what we had in the death
# simulation: each individual had the same probability, a Bernoulli trial
# occurred for each individual, and the number dead is the sum across
# individuals. In other words, when we scaled up from the individual-level
# stochastic process, we got a new but known stochastic process at the larger
# scale.

# Functions for distributions in Python are in the SciPy library

from scipy.stats import binom

# Now we can plot the theoretical PMF for the binomial against the outcome of
# the individual-based simulations. They are very close.

X = np.arange(0, 45)
pmf = binom.pmf(X, N, p_d)
plt.plot(X, pmf, marker='o', linestyle='none', fillstyle='none')
plt.title("Binomial distribution (points) vs simulation (bars)")
plt.show()

# The binomial emerges from the underlying (biological) stochastic process!

# Mean (expected value)

# Theoretical E(X) = Np
print(N * p_d)
# Observed
print(np.mean(N_dead))

# Variance

# Theoretical Var[X] = Np(1-p)
print(N * p_d * (1 - p_d))
# Observed (sample variance)
print(np.var(N_dead, ddof=1))
