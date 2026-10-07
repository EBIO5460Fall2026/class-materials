# This was our Poisson Ricker model from last week

Ricker_Poisson_base <- function(Nt, beta, m, alpha){
    # Births
    births <- rpois( length(Nt), Nt * beta )
    # Density dependent and density independent survival
    survivors <- rbinom( length(Nt), births, (1 - m) * exp( -1 * alpha * Nt ) )
    return(survivors)
}

# We need to make some tweaks to these models for confronting with data. It's
# not possible to distinguish between beta and m from data. This is called
# non-identifiability. We can change one or the other for identical effect, not
# just in the expectation (mean of the stochastic model) but also in the full
# distribution. This is a consequence of how the different stochastic processes
# combine together in this case. Beta and m form a composite parameter, which is
# the finite, density-independent growth rate of the population:
#
#     R = beta ( 1 - m )


# Poisson births, Pois(beta), combined with binomial, density-independent
# survival, Binom(births, 1-m), gives Poisson survivors. We'll next check that
# this is indeed true.

# Vectorized simulation

sims <- 1000000
Nt <- rep(25, sims)
beta <- 5
m <- 0.2
births <- rpois(length(Nt), Nt * beta)
di_survivors <- rbinom(length(Nt), births, 1 - m)

# Distribution (PMF) of di_survivors
hist(di_survivors, breaks=seq(min(di_survivors)-0.5, max(di_survivors)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="number surviving",
    ylab="Probability mass", main="")

# Show that this two-step birth survival is Poisson with mean N_t R

R = beta * (1 - m)

# Add theoretical PMF to plot
X <- min(di_survivors):max(di_survivors)
pmf <- dpois(X, Nt[1] * R)
points(X, pmf, col="darkorange")
title("Poisson(Nt R) (points) vs simulation (bars)")

# Next we account for density-dependent survival. This is less biologically
# intuitive because density-independent survival and density-dependent survival
# actually happen at the same time, not one after the other. The model, while
# identical in outcome, is not following the event sequencing of the biological
# algorithm. But this allows us to combine the non-identifiable parameters.

alpha <- 0.002
survivors <- rbinom( length(Nt), di_survivors, exp( -1 * alpha * Nt ) )

# Distribution (PMF) of survivors
hist(survivors, breaks=seq(min(survivors)-0.5, max(survivors)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="number surviving",
    ylab="Probability mass", main="")

# Show that this overall is Poisson with mean N_t R exp(-alpha * Nt)

mu <- Nt[1] * R * exp(-alpha * Nt[1])

# Add theoretical PMF to plot
X <- min(survivors):max(survivors)
pmf <- dpois(X, mu)
points(X, pmf, col="darkorange")
title("Poisson(Nt R exp(-alpha Nt)) (points) vs simulation (bars)")


# Combining this into a (much neater) function

Ricker_Poisson <- function(Nt, R, alpha){

    # Births and density independent survival
    births_di_survive <- rpois( length(Nt), Nt * R )

    # Density dependent survival
    survivors <- rbinom( length(Nt), births_di_survive, exp( -1 * alpha * Nt ) )

    return(survivors)
}

# Vectorized simulation
sims <- 1000000
Nt <- rep(25, sims)
R <- 5 * ( 1 - 0.2 )
alpha <- 0.002
Ntp1 <- Ricker_Poisson(Nt, R, alpha)

# Distribution (PMF) of Ntp1
hist(Ntp1, breaks=seq(min(Ntp1)-0.5, max(Ntp1)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="number surviving",
    ylab="Probability mass", main="")

# Show that this overall is Poisson with mean N_t R exp(-alpha * Nt)

mu <- Nt[1] * R * exp(-alpha * Nt[1])

# Add theoretical PMF to plot
X <- min(survivors):max(survivors)
pmf <- dpois(X, mu)
points(X, pmf, col="darkorange")
title("Poisson(Nt R exp(-alpha Nt)) (points) vs simulation (bars)")


# Deriving the Negative-binomial-binomial-gamma Ricker
# This is the most complex of the models with all sources of stochasticity

# This function takes an input vector N_{t}, abundance at time t, and outputs
# N_{t+1}, abundance at time t+1. We add the following stochastic parameters:
#   kD: heterogeneity in individual R
#   kE: environmental stochasticiy in R
#   p:  probability of female offspring
Ricker_NBBG <- function(Nt, R, alpha, kD, kE, p=0.5) {

    # Environmental stochasticity: heterogeneity in R (birth rate/DI mortality)
    # between times or locations (gamma with mean R)
    Rtx <- rgamma(length(Nt), shape=kE, scale=R/kE)

    # Heterogeneity in sex (binomial)
    females <- rbinom(length(Nt), Nt, p)

    # Heterogeneity in individual birth rate (gamma) plus density-independent
    # survival (R = births (1-mortality) (1/p); mortality is binomial, so
    # compound distribution is negative binomial with mean (1/p) R Nt).
    # (add 1 to size when females=0 to avoid NaN's)
    births <- rnbinom( length(Nt), size = kD * females + (females==0),
            mu = (1/p) * females * Rtx )

    # Density dependent survival
    survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )

    return(survivors)
}


# Vectorized simulation of this stochastic model
sims <- 1000000
Nt <- rep(25, sims)
R <- 5 * ( 1 - 0.2 )
alpha <- 0.002
kD <- 10
kE <- 10
p <- 0.5
Ntp1 <- Ricker_NBBG(Nt, R, alpha, kD, kE, p)

# Distribution (PMF) of N_{t+1}
hist(Ntp1, breaks=seq(min(Ntp1)-0.5, max(Ntp1)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="number surviving",
    ylab="Probability mass", main="")


# Show that the theoretical PMF is correct

# The function for the PMF is quite a beast with summations and integrals
# calculated numerically. We won't consider the technical details. It's in the
# source directory. The point is, despite being complicated, the PMF can be
# calculated accurately.
source("source/Ricker_dnbinombinomgamma.R")

# Add theoretical PMF to plot
X <- min(Ntp1):max(Ntp1)
pmf <- Ricker_dnbinombinomgamma(X, rep(Nt, length(X)), p, R, alpha, kD, kE)
points(X, pmf, col="darkorange")
title("NBBG(Nt R exp(-alpha Nt)) (points) vs simulation (bars)")
