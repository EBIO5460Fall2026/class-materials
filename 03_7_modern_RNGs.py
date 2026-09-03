# Modern random number generators (RNGs) in Python

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Now that we've seen a full cycle of hand-coded algorithms through to
# agent-based simulations, this is a good point to discard our hand-coded RNG
# and replace it with modern, ultra-reliable tools from the numpy library. One
# can also do random generation with the SciPy library but Numpy has more modern
# random number generation.

# Numpy Uniform(0,1) random numbers

# Initiate the random number generator

rng = np.random.default_rng(98765)

# The object `rng` will keep track of its state. This is a slightly different
# mechanism for keeping track from what I illustrated in our hand-coded RNG but
# the concept is the same. The object `rng` is set up in such a way that keeping
# track is an internal process. It is using an object-oriented programming
# approach. Essentially `rng` "knows" what its state is as we continue to use
# it.

# Now use it

U = rng.random(10)
print(U)

# Use it again

U = rng.random(10)
print(U)             #different set


# Make a stand-alone reproducible code chunk by resetting the seed

rng = np.random.default_rng(736)

draws = 10000

X = rng.random(size=draws)
print(X)
print(np.mean(X)) #Expected = 0.5

plt.figure()
plt.hist(X, bins=30, density=True)
plt.xlabel("x")
plt.ylabel("Density")
plt.plot([0, 1], [1, 1])
plt.title("Uniform(0,1) distribution (line) vs random draws (histogram)")
plt.show()


# Numpy Bernoulli random numbers (via Binomial with number of trials = 1)

# Initiate the random number generator (optional if previously initiated)

rng = np.random.default_rng(631)

draws = 10000
p = 0.3

X = rng.binomial(n=1, p=p, size=draws)

# Notice argument names: n is the number of trials in the binomial distribution,
# whereas size is the number of draws. Python and R have these names reversed,
# so watch out!

print(X)
print(np.mean(X)) #Expected = p

plt.figure()
plt.hist(X, bins=[-0.5, 0.5, 1.5], density=True)
plt.xticks([0, 1])
plt.xlabel("x")
plt.ylabel("Relative frequency")
plt.plot([0, 1], [1-p, p], marker='o', linestyle='none', fillstyle='none')
plt.title("Bernoulli distribution (points) vs random draws (histogram)")
plt.show()


# Numpy Binomial random numbers

# Initiate the random number generator
rng = np.random.default_rng(2293)

draws = 100000
n = 100    #number of Bernoulli trials
p = 0.2

X = rng.binomial(n, p, size=draws)
print(X)
print(np.mean(X)) #Expected = np

bins = np.arange(np.min(X) - 0.5, np.max(X) + 0.5 + 1)
x = np.arange(np.min(X), np.max(X) + 1)
pmf = stats.binom.pmf(x, n, p)
plt.figure()
plt.hist(X, bins=bins, density=True)
plt.plot(x, pmf, marker='o', linestyle='none', fillstyle='none')
plt.xlabel("x")
plt.ylabel("Relative frequency")
plt.title("Binomial distribution (points) vs random draws (histogram)")
plt.show()
