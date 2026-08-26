# Generating pseudo-random numbers
# We model a Uniform(0,1) Data Generating Process
# Classic Lewis et al. 1969 algorithm

# This algorithm can sometimes experience integer overflow because R's default
# integer type is 32-bit but it's usually fine. Demonstration only.

n <- 10  # how many random numbers?

a <- 16807
m <- 2147483647

# Random integers

I <- 6
for (j in 1:n) {
  I <- (a * I) %% m
  print(I)
}

# Random "continuous" numbers in 0-1

I <- 6
for (j in 1:n) {
  I <- (a * I) %% m
  U <- I / m
  print(U)
}

# Turn this into a function returning n uniform pseudo-random numbers

rand_unif <- function(n, seed) {
  a <- 16807
  m <- 2147483647

  U <- numeric(n)

  I <- seed
  for (j in 1:n) {
    I <- (a * I) %% m
    U[j] <- I / m
  }

  U
}

# Use the function to generate a large vector of random numbers

U <- rand_unif(n = 1000000, seed = 6)

print(U)

# Histogram

hist(
  U,
  breaks = 20,
  xlab = "U",
  ylab = "Frequency",
  main = ""
)

# In stochastic simulation, the "histogram is the distribution"
# Histogram as density (proper distribution)
# Density is frequency divided by the total area under the curve
# This normalizes frequency to a density because now the total area
# under the curve = 1

hist(
  U,
  breaks = 20,
  probability = TRUE,
  xlab = "U",
  ylab = "Density",
  main = ""
)

abline(h = 1.0, col = "red", lty = 2)

# As n -> Inf, the simulation approaches "truth"

# In this algorithm, the seed cannot be zero (or the values are all zero)
# and the sequence never visits zero

rand_unif(n = 10, seed = 0)
