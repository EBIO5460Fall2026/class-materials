# Continuous-time birth-death IBM

# Simulation control
set.seed(95756)
t <- 0
t_max <- 30

# Initial biological state

# --individual scale
N_0 <- 10
alive <- rep(1, N_0)
b <- rep(0.25, N_0)    # birth rate
m <- rep(0.2, N_0)     # death rate

# --population scale
N_t <- sum(alive)      # scale transition (sum over individuals -> population)
N <- N_t               # track N
times <- t             # track time

# Simulate events
while ( t < t_max & N_t > 0 ) {

    # calculate event intensity
    total_birth <- sum(b * alive)
    total_death <- sum(m * alive)
    intensity <- total_birth + total_death

    # waiting time to next event
    w <- rexp(1, rate=intensity)

    # advance time
    t <- t + w

    # which event type occurs?
    event_birth <- rbinom(1, size=1, prob=total_birth / intensity)

    # which individual does the event happen to?
    i <- sample(which(alive == 1), size=1)  #equal p, identical individuals

    # update system state (do event)
    if ( event_birth == 1 ) {

        # add a new individual and its traits
        alive <- c(alive, 1)
        b <- c(b, 0.25)
        m <- c(m, 0.2)

    } else {

        # death
        alive[i] <- 0

    }

    N_t <- sum(alive)
    N <- c(N, N_t)
    times <- c(times, t)
}

print(alive)

plot(times, N, type="s", ylim=c(0, max(N)), col="blue", xlab="Time", ylab="N",
    main="Birth-death IBM")
