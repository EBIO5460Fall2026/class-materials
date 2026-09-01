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

# 00 individuals
# Probability of death in a year = 0.2
# How many die in the first year?

N = 100
alive = np.full(N, 1)    #alive = 1
p_death = np.full(N, 0.2)

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

    # Simulate the death process for each individual

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
    # if i % 1000 == 0:  # monitoring
    #     print(i)

print(alive)   #each simulation is one column

# Each individual has its own Bernoulli distribution
# Extract all the outcomes for individual 37

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

# See, it emerges from the process!
# Now go see the slide for Binomial emerging


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





# Modern random number generators (RNGs)

# Now that we've seen a full cycle of hand-coded algorithms through to
# agent-based simulations, this is a good point to discard our hand-coded RNG
# and replace it with modern, ultra-reliable tools from the numpy library. One
# can also do random generation with the SciPy library but Numpy has more modern
# random number generation.

# Numpy Uniform(0,1) random numbers

# Initiate the random number generator
# The object `rng` will keep track of its state

rng = np.random.default_rng(98765)

# Now use it

U = rng.random(10)
print(U)

# Use it again

U = rng.random(10)
print(U)             #different set


# Numpy Bernoulli random numbers (via Binomial with trials=1)

# Initiate the random number generator (optional if previously initiated)

rng = np.random.default_rng(631)

draws = 10000
p = 0.3

X = rng.binomial(1, p, size=draws)
print(X)
print(np.mean(X)) #Expected = p


# Numpy Binomial random numbers

# Initiate the random number generator
rng = np.random.default_rng(2293)

draws = 100000
n = 100    #number of Bernoulli trials
p = 0.2

X = rng.binomial(n, p, size=draws)
print(X)
print(np.mean(X)) #Expected = np


#--------------------------------------------------------------------------------
# Exponential distribution and event waiting times
#--------------------------------------------------------------------------------

c = 0.02         #roughly probability of death per year (so, lives ca 50 years)
delta_t = 1/365  #years
p = c * delta_t  #probability per day
print(p)

state = 0
t = 0
while state == 0:
    t = t + delta_t
    if ( rng.binomial(1, p, 1) ):
        state = 1
print(t)


# Turn this algorithm into a function

def death_continous_time(c, delta_t):
    p = c * delta_t
    alive = 1
    t = 0
    while alive:
        t = t + delta_t
        if ( rng.binomial(1, p, 1) ):
            alive = 0
    return t


# Conduct a simulation study

rng = np.random.default_rng(58773)

c = 0.02
delta_t = 1/365
reps = 1000
time_to_death = np.empty(reps)

for i in range(reps):
    time_to_death[i] = death_continous_time(c, delta_t)

print(time_to_death)
print(np.mean(time_to_death)) #Expected: 1 / c
print(np.var(time_to_death, ddof=1))  #well, that's unrealistic

plt.figure()
plt.hist(time_to_death, bins=20, density=True)
plt.xlabel("Time to death")
plt.ylabel("Density")
plt.show()

# Expected
t = np.linspace(0, np.max(time_to_death) + 1, 100)
pdf = c * np.exp(-c * t)
plt.plot(t, pdf)

