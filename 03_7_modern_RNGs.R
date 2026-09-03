# Modern random number generators (RNGs) in R

# Now that we've seen a full cycle of hand-coded algorithms through to
# agent-based simulations, this is a good point to discard our hand-coded RNG
# and replace it with modern, ultra-reliable tools. In R, these are built in to
# the base library.

# Uniform(0,1) random numbers

# Initiate the random number generator

set.seed(98765)

# The RNG in R keeps track of its state as a global object. The state of the RNG
# is updated automatically.

# Now use it

U <- runif(10)
U

# Use it again

U <- runif(10)
U             #different set


# Make a stand-alone reproducible code chunk by resetting the seed

set.seed(736)

draws <- 10000

X <- runif(draws)
head(X)
mean(X) #Expected = 0.5

hist(X, breaks=30, probability=TRUE, xlab="x", ylab="Density", col="lightblue",
    main="Uniform(0,1) distribution (line) vs random draws (histogram)")
lines(c(0, 1), c(1, 1), col="darkorange")


# Bernoulli random numbers (via Binomial with number of trials = 1)

# Initiate the random number generator (optional if previously initiated)

set.seed(631)

draws <- 10000
p <- 0.3

X <- rbinom(n=draws, size=1, prob=p)

# Notice argument names: n is the number of draws, whereas size is the number of
# trials in the binomial distribution. Python and R have these names
# reversed, so watch out!

head(X)
mean(X) #Expected = p

hist(X, breaks=c(-0.5, 0.5, 1.5), probability=TRUE, xaxt="n", xlab="x",
    col="lightblue", ylab="Relative frequency",
    main="Bernoulli distribution (points) vs random draws (histogram)")
axis(1, at=c(0, 1))
points(c(0, 1), c(1 - p, p), col="darkorange")


# Binomial random numbers

# Initiate the random number generator

set.seed(2293)

draws <- 100000
n <- 100    #number of Bernoulli trials
p <- 0.2

X <- rbinom(draws, size=n, prob=p)

head(X)
mean(X) #Expected = np

bins <- seq(min(X) - 0.5, max(X) + 0.5, by=1)
x <- min(X):max(X)
pmf <- dbinom(x, size=n, prob=p)

hist(X, breaks=bins, probability=TRUE, xlab="x", ylab="Relative frequency", col="lightblue",
    main="Binomial distribution (points) vs random draws (histogram)")
points(x, pmf, col="darkorange")
