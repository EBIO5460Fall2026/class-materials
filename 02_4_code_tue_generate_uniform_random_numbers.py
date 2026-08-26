# Generating pseudo-random numbers
# We model a Uniform(0,1) Data Generating Process
# Classic Lewis et al 1969 algorithm

n = 10 #how many random numbers?

a = 16807
m = 2147483647

# Random integers

I = 6
for j in range(n):
    I = (a * I) % m
    print(I)

# Random "continuous" numbers in 0-1

I = 6
for j in range(n):
    I = (a * I) % m
    U = I / m
    print(U)


# Turn this into a function returning n uniform pseudo-random numbers

import numpy as np

def rand_unif(n, seed):
    a=16807
    m=2147483647
    U = np.empty(n)

    I = seed
    for j in range(n):
        I = (a * I) % m
        U[j] = I / m

    return U

# Use the function to generate a large vector of random numbers

U = rand_unif(n=1000000, seed=6)
print(U)

# Histogram

import matplotlib.pyplot as plt

plt.figure()
plt.hist(U, bins=20)
plt.xlabel("U")
plt.ylabel("Frequency")
plt.show()

# In stochastic simulation, the "histogram is the distribution"
# Histogram as density (proper distribution)
# Density is frequency divided by the total "area under the curve"
# This normalizes frequency to a density because now the total area under the
# curve = 1

plt.figure()
plt.hist(U, bins=20, density=True)
plt.axhline(1.0, color="red", linestyle="--")
plt.xlabel("U")
plt.ylabel("Density")
plt.show()

# As n -> Inf, the simulation approaches "truth" (red line)

# In this algorithm, the seed cannot be zero (or the values are all zero) and
# the sequence never visits zero

rand_unif(n=10, seed=0)
