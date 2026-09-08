import numpy as np
import matplotlib.pyplot as plt

### RNG from last time

def rand_unif(n, seed):
    a = 16807
    m = 2147483647
    U = np.empty(n)

    I = seed
    for j in range(n):
        I = (a * I) % m
        U[j] = I / m

    return U, I


### Start the RNG

U, I_state = rand_unif(1, seed=74732)


### Individual-based model

# 100 individuals
# Probability of death in a year = 0.2
# How many die in the first year?

N = 100
alive = np.full(N, 1)    #alive = 1
p_death = np.full(N, 0.2)

# Simulate death for each individual
for i in range(N):
    U, I_state = rand_unif(1, I_state)
    if U <= p_death[i]:
        alive[i] = 0    #individual dies

print(sum(alive))
print(sum(np.logical_not(alive))) #number dead


### Package into a function

def death_IBM(N, p_d, I_state):

    alive = np.full(N, 1)        #alive = 1
    p_death = np.full(N, p_d)

    for i in range(N):
        U, I_state = rand_unif(1, I_state)
        if U <= p_death[i]:
            alive[i] = 0     #individual dies

    return alive, I_state

# Check that it works

alive, I_state = death_IBM(N=100, p_d=0.2, I_state=I_state)
print(alive)
print(sum(alive))


### Simulate IBM repeatedly

U, I_state = rand_unif(1, seed=924505)

# Simulation control
sims = 100000

# Ecological parameters
N = 100
p_d = 0.2

# Storage for results
alive = np.empty((N, sims), dtype=int)

# Computation
for i in range(sims):
    alive[:,i], I_state = death_IBM(N, p_d, I_state)
    #if i % 1000 == 0:  # uncomment for monitoring
    #    print(i)

# Output
print(alive)   #each simulation is one column


### Bernoulli distribution

# Extract all realizations for one individual
individual = 37
dead = (alive[individual,:] - 1) * -1 #convert alive to dead

# Plot histogram
plt.figure()
plt.hist(dead, bins=[-0.5, 0.5, 1.5], density=True)
plt.xticks([0, 1])
plt.xlabel("Dead")
plt.ylabel("Relative frequency")

# Compare to theoretical PMF
plt.plot([0, 1], [1-p_d, p_d], marker='o', linestyle='none', fillstyle='none')
plt.title("Bernoulli distribution (points) vs simulation (bars)")

### Population-level outcome

# Total alive or dead (sum columns)
N_alive = np.sum(alive, axis=0)
N_dead = np.sum(np.logical_not(alive), axis=0)
print(N_alive)
print(N_dead)

# Distribution (PMF) of number dead
plt.figure()
bins = np.arange(np.min(N_dead)-0.5, np.max(N_dead)+0.5+1)
plt.hist(N_dead, bins, density=True)
plt.xlabel("Number dead")
plt.ylabel("Probability mass")


### Binomial distribution

from scipy import stats

# Add theoretical PMF to plot
X = np.arange(np.min(N_dead), np.max(N_dead)+1)
pmf = stats.binom.pmf(X, N, p_d)
plt.plot(X, pmf, marker='o', linestyle='none', fillstyle='none')
plt.title("Binomial distribution (points) vs simulation (bars)")

# Theoretical mean (expected value):
print(N * p_d)

# Simulation mean:
print(np.mean(N_dead))

# Theoretical variance:
print(N * p_d * (1 - p_d))

# Simulation sample variance:
print(np.var(N_dead, ddof=1))
