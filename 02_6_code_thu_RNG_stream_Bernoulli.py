import numpy as np
import matplotlib.pyplot as plt

# How a random number generator stream works:
# the current state of the generator is kept to start the next batch.
# To illustrate, we can modify our previous function to also return the state
# of the generator (add `I` to the `return` statement).

def rand_unif(n, seed):
    a = 16807
    m = 2147483647
    U = np.empty(n)

    I = seed
    for j in range(n):
        I = (a * I) % m
        U[j] = I / m

    return U, I          #also return generator state I


# Typically, a random number generator is initialized at first

U, I_state = rand_unif(n=1, seed=9033)

# Now we can continue to generate from the same stream of numbers without
# supplying a new fixed seed, by feeding back in the current state

U, I_state = rand_unif(10, I_state)
print(U)

# With most current inbuilt RNGs we won't need to do this manually but this
# illustrates conceptually what is happening with an RNG stream.

# PS. Also, you won't ever want to use the RNG I'm using here. There are better
# (satisfying more stringent tests of randomness), faster, and more convenient
# ones in the NumPy library. But it's not merely a toy either. This algorithm
# was used widely in scientific computing until the early 2000s. I'm using it
# only for demonstration of algorithm concepts. This and the remainder of the
# usage below is to give an explicit algorithmic sense of what is going on in
# a stochastic simulation on the computer. 


# Binary stochastic events

# Conceptualize a biological process, such as birth, death, disturbance,
# mutation etc, as an event that occurs with probability p. There are only two
# possible outcomes, e.g. the event occurs or doesn't occur, or alternatively
# the event occurs or some other event occurs.

# How can we generate a binary stochastic event with probability p?

# Algorithm
# generate a uniform random number, U
# if U <= p
#     the event occurs
# else
#     it doesn't

p = 0.3

U, I_state = rand_unif(1, I_state)

if U <= p:
    event = True
else:
    event = False

# Event occurred?

print(U)
print(event)


# Bernoulli distribution

# The above stochastic process gives rise to what is known as the Bernoulli
# distribution. We can simulate multiple binary stochastic events to recover the
# PMF and other properties (mean, variance) of the Bernoulli distribution.

# Repeated realizations of a binary stochastic event
# Vectorized

p = 0.3
reps = 10000

U, I_state = rand_unif(reps, I_state)

# True if event occurs, False otherwise

events = U <= p
print(events)

# Number of occurrences (numpy considers True = 1, False = 0)

count = sum(events)
print(count)

# Mean (expected value)

# Theoretical E(X) = p
print(p)
# Observed
print(count / reps)
# or
print(np.mean(events))

# Variance

# Theoretical Var[X] = p(1-p)
print(p * (1 - p))
# Observed (sample variance)
print(np.var(events, ddof=1))

# Frequency histogram

# matplotlib doesn't automatically treat False=0, True=1
events_int = events.astype(int)
plt.figure()
# bins centered on integer with width = 1
plt.hist(events_int, bins=[-0.5, 0.5, 1.5])
plt.xticks([0, 1])
plt.xlabel("Event")
plt.ylabel("Frequency")
plt.show()

# Histogram of the relative frequency to represent probability mass

# A probability mass is really a point mass, so conceptually there are not bins
# with a width spanning a range of values (e.g. since the only possible values
# are the integers 0 and 1, the probability that X = 1.01 is zero). Hence a
# density histogram (density=TRUE) is not the technically correct visualization.
# Many people represent this as a bar plot with thin bars or a lollipop plot
# with lines and a point for the point mass. But it's also common to use a
# density histogram function with bin width equal to 1, which gives the correct
# relative frequencies to compare the estimated probability mass of each
# integer. This is my preferred plot because the adjacent bars are easier to
# compare in relative length, area is a strong channel to represent mass, it's
# less code, and the code is consistent with code for continuous distributions.
# However, it's important to remember the distinction between density for
# continuous distributions and point mass for discrete distributions.

plt.figure()
plt.hist(events_int, bins=[-0.5, 0.5, 1.5], density=True)
plt.xticks([0, 1])
plt.xlabel("Event")
plt.ylabel("Relative frequency")
plt.show()

# Add the theoretically expected values from the Bernoulli distribution to see
# that the simulation does indeed match.

plt.plot([0, 1], [1-p, p], marker='o', linestyle='none', fillstyle='none')
plt.title("Bernoulli distribution (points) vs simulation (bars)")
plt.show()

# Agent-based and individual-based simulation

# We now have all the tools we need to do our first agent-based simulation!
# We'll simulate a simple death process.

# Imagine a population with 100 individuals (e.g. 100 Tribolium beetles in my
# lab, 100 plants, 100 parasites, 100 lemurs)
# Each has probability of death in a year equal to 0.2
# How many die in the first year?

# Start the RNG with a new seed to make the simulation reproducible

U, I_state = rand_unif(1, seed=74732)

# In an agent-based simulation, we want to keep track of each agent separately.
# Here, our agents are individuals. So we set up data structures to keep track
# of all the attributes of each individual and we initiate the simulation by
# giving each individual its attributes. For this simple simulation, the only
# attributes are each individual's "alive" status and probability of death.

N = 100
alive = np.full(N, 1)    #alive = 1
p_death = np.full(N, 0.2)

# In this simple scenario, we could vectorize the simulation (recommended for
# speed, simplicity, and code readability) but to be general and allow extension
# to more complicated scenarios where vectorization is impossible, we'll iterate
# explicitly over each individual.

# Simulate the death process for each individual

for i in range(N):
    U, I_state = rand_unif(1, I_state)
    if U <= p_death[i]:
        alive[i] = 0    #individual dies

print(alive)    #which individuals are dead and alive
print(sum(alive))    #85 remaining alive
print(sum(np.logical_not(alive)))    #15 dead

# We can run that chunk of code (minus setting the seed) multiple times to
# simulate different realizations of the stochastic process. In each run,
# different individuals die. Furthermore, different numbers of individuals die.
# There is a new stochastic process that emerges at the scale of the population
# for the number of individuals that die. In other words, when we scale up from
# an individual or unit-level stochastic process, we get a new stochastic
# process at the larger scale.
