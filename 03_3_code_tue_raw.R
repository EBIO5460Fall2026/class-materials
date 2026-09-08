
### RNG from last time

rand_unif <- function(n, seed) {
    a <- 16807
    m <- 2147483647
    U <- numeric(n)

    I <- seed
    for ( j in 1:n ) {
        I <- (a * I) %% m
        U[j] <- I / m
    }

    list(U=U, I_state=I)
}


### Start the RNG

rng_out <- rand_unif(1, seed = 74732)


### Individual-based model

# 100 individuals
# Probability of death in a year = 0.2
# How many die in the first year?

N <- 100
alive <- rep(1, N)    #alive = 1
p_death <- rep(0.2, N)

# Simulate death for each individual
for ( i in 1:N ) {
    rng_out <- rand_unif(1, rng_out$I_state)
    U <- rng_out$U
    if ( U <= p_death[i] ) {
        alive[i] <- 0    #individual dies
    }
}

alive
sum(!alive) #number dead


### Package into a function

death_IBM <- function(N, p_d, I_state) {

    alive <- rep(1, N)        #alive = 1
    p_death <- rep(p_d, N)

    for ( i in 1:N ) {
        rng_out <- rand_unif(1, I_state)
        U <- rng_out$U
        I_state <- rng_out$I_state
        if ( U <= p_death[i] ) {
            alive[i] <- 0    #individual dies
        }
    }

    return(list(alive=alive, I_state=I_state))
}

# Check that it works

sim_out <- death_IBM(N=100, p_d=0.2, I_state=rng_out$I_state)
alive <- sim_out$alive
alive
sum(alive)


### Simulate IBM repeatedly

rng_out <- rand_unif(1, seed=924505)
I_state <- rng_out$I_state

# Simulation control
sims <- 100000

# Ecological parameters
N <- 100
p_d <- 0.2

# Storage for results
alive <- matrix(NA, nrow=N, ncol=sims)

# Computation
for ( i in 1:sims ) {
    sim_out <- death_IBM(N, p_d, I_state)
    alive[,i] <- sim_out$alive
    I_state <- sim_out$I_state
    #if ( i %% 1000 == 0 ) {  # uncomment for monitoring
    #    print(i)
    #}
}

# Output
alive[1:10, 1:10]   #each simulation is one column


### Bernoulli distribution

# Extract all realizations for one individual
individual <- 38
dead <- (alive[individual,] - 1) * -1 #convert alive to dead

# Plot histogram
hist(dead, breaks=c(-0.5, 0.5, 1.5), probability=TRUE, col="lightblue",
    axes=FALSE, ann=FALSE)
axis(1, at=c(0,1), lwd=0)
axis(2)
title(xlab="Dead", ylab="Relative frequency")

# Compare to theoretical PMF
points(c(0,1), c(1-p_d, p_d), col="darkorange")
title(main="Bernoulli distribution (points) vs simulation (bars)")


### Population-level outcome

# Total alive or dead (sum columns)
N_alive <- colSums(alive)
N_dead <- colSums(!alive)
head(N_alive)
head(N_dead)

# Distribution (PMF) of number dead
hist(N_dead, breaks=seq(min(N_dead)-0.5, max(N_dead)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="Number dead",
    ylab="Probability mass", main="")


### Binomial distribution

# Add theoretical PMF to plot
X <- min(N_dead):max(N_dead)
pmf <- dbinom(X, N, p_d)
points(X, pmf, col="darkorange")
title("Binomial distribution (points) vs simulation (bars)")

# Theoretical mean (expected value):
N * p_d

# Simulation mean:
mean(N_dead)

# Theoretical variance:
N * p_d * (1 - p_d)

# Simulation sample variance:
var(N_dead)
