# Stochastic Ricker model
# Assumptions: fish lay eggs at the start, adults cannabilize eggs

# Individual perspective

# IBM algorithm

# Individual i gives birth to eggs (Poisson)

birth_ind_i <- function(beta) {
    B <- rpois(1, beta)
    return(B)
}

eggs_i <- birth_ind_i(5)

# Eggs survive cannibalism (Bernoulli) and intrinsic mortality (Bernoulli)

# Each egg is eaten with probability 1 - exp(-alpha* Nt) or dies with
# probability m

egg_surv <- function(eggs, m, alpha, N) {

    surviving_eggs <- 0
    for ( egg in 1:eggs ) {
        survive <- TRUE
        # cannibalized?
        if ( rbinom(1, 1, (1 - exp(-alpha * N)) ) ) {
            survive <- FALSE
        }
        # intrinsic death?
        if ( rbinom(1, 1, m ) ) {
            survive <- FALSE
        }
        # record surviving egg
        if ( survive ) {
            surviving_eggs <- surviving_eggs + 1
        }
    }
    return(surviving_eggs)
}

# Scaling up from individual to population (sum over individuals)

Ricker_IBM <- function(Nt, beta, m, alpha) {

    Ntp1 <- 0
    # for each adult
    for ( i in seq_len(Nt) ) {
        eggs_i <- birth_ind_i(beta)
        surv_offspring <- egg_surv(eggs_i, m, alpha, Nt)
        Ntp1 <- Ntp1 + surv_offspring
    }
    return(Ntp1)
}

Nt <- 20
beta <- 5
m <- 0.2
alpha <- 0.002

Ntp1 <- Ricker_IBM(Nt, beta, m, alpha)


# Simulation study 1 (full IBM)

set.seed(55157)

sims <- 10000
Nt <- 20
beta <- 5
m <- 0.2
alpha <- 0.002

Ntp1 <- rep(NA, sims)
for ( i in 1:sims) {
    Ntp1[i] <- Ricker_IBM(Nt, beta, m, alpha)
}

# Distribution (PMF) of Ntp1
hist(Ntp1, breaks=seq(min(Ntp1)-0.5, max(Ntp1)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="N_{t+1}",
    ylab="Probability mass", main="")

# Add theoretical PMF (Poisson) to plot
X <- min(Ntp1):max(Ntp1)
pmf <- dpois(X, Nt * beta * (1 - m) * exp(-alpha * Nt))
points(X, pmf, col="darkorange")
title("Poisson distribution (points) vs simulation (bars)")


# Simulation study 2 (equivalent vectorized model)

# Simplified vectorized version

Ricker_Poisson <- function(Nt, beta, m, alpha){
    # Births
    births <- rpois( length(Nt), Nt * beta )
    # Density dependent and density independent survival
    survivors <- rbinom( length(Nt), births, (1 - m) * exp( -1 * alpha * Nt ) )
    return(survivors)
}


set.seed(3098)

sims <- 100000
Nt <- 20
beta <- 5
m <- 0.2
alpha <- 0.002

Ntp1 <- rep(NA, sims)
for ( i in 1:sims) {
    Ntp1[i] <- Ricker_Poisson(Nt, beta, m, alpha)
}


# Distribution (PMF) of Ntp1
hist(Ntp1, breaks=seq(min(Ntp1)-0.5, max(Ntp1)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="N_{t+1}",
    ylab="Probability mass", main="")

# Add theoretical PMF to plot
X <- min(Ntp1):max(Ntp1)
pmf <- dpois(X, Nt * beta * (1 - m) * exp(-alpha * Nt))
points(X, pmf, col="darkorange")
title("Poisson distribution (points) vs simulation (bars)")


#Simulate one time series
N_0 <- 20
beta <- 5
m <- 0.2
alpha <- 0.002

t <- 0:200
N <- t * NA
N[1] <- N_0
for (i in 1:max(t)) {
    N[i+1] <- Ricker_Poisson(N[i], beta, m, alpha)
}
plot(t, N, type="l")
