# NBBG model with jags

library(R2jags)
library(bayesplot)
library(ggplot2)

tribdata <- read.csv("data/ricker_data.csv")
Nt <- round(tribdata$At)
Ntp1 <- round(tribdata$Atp1)

n = length(Nt)
data <- list(
    Nt = Nt,
    Ntp1 = Ntp1,
    n = n
)

# We need good starting values for this model, especially for the hidden nodes,
# births and females
inits <- function() {
    list(
        lnR = log(1),
        lnalpha = log(0.001),
        lnkD = log(10),
        lnkE = log(10),
        births = Ntp1
       # females = pmax(1, round(Nt/2))
    )
}

params <- c("R", "alpha", "kD", "kE")

# Initial simulation to check chain behavior
Ricker_nbbg_sim <- jags(
    data = data,
    inits = inits,
    parameters.to.save = params,
    model.file = "ricker_nbbg.jags",
    n.chains = 1,
    n.iter = 4000,
    n.burnin = 0,      #keep all
    n.thin = 1         #iterations
)
Ricker_nbbg_sim

samples <- as.data.frame(Ricker_nbbg_sim$BUGSoutput$sims.matrix)
mcmc_trace(samples, "R")
mcmc_trace(samples, "R", window=c(0,500))
mcmc_trace(samples, "alpha")
mcmc_trace(samples, "alpha", window=c(0,500))
mcmc_trace(samples, "kD")
mcmc_trace(samples, "kD", window=c(0,500))
mcmc_trace(samples, "kE")
mcmc_trace(samples, "kE", window=c(0,500))
mcmc_acf(samples, c("R","alpha","kD","kE"))


# Multiple chains now and bumping up the iterations
Ricker_nbbg_sim <- jags(
    data = data,
    inits = inits,
    parameters.to.save = params,
    model.file = "ricker_nbbg.jags",
    n.chains = 4,
    n.iter = 40000,
    n.burnin = 30000,
    n.thin = 2
)
# Rhat is struggling to show great convergence and n.eff is very low but this
# could be acceptable with more computation
Ricker_nbbg_sim

# Check convergence
# The chains on the other hand show trouble with mixing
samples_array <- Ricker_nbbg_sim$BUGSoutput$sims.array
mcmc_trace(samples_array, "R")
mcmc_trace(samples_array, "R", window=c(0,500))
mcmc_trace(samples_array, "alpha", window=c(0,500))
mcmc_trace(samples_array, "kD", window=c(0,500))
mcmc_trace(samples_array, "kE", window=c(0,500))

# alpha is mixing poorly
mcmc_trace(samples_array, "alpha", window=c(4500,5000))
mcmc_trace(samples_array, "R", window=c(4500,5000))

# Posterior distributions (including 2D view)
samples <- as.data.frame(Ricker_nbbg_sim$BUGSoutput$sims.matrix)
head(samples)
mcmc_hist(samples, pars=c("R","alpha","kD","kE"))
mcmc_pairs(samples, pars=c("R","alpha","kD","kE"))
mcmc_pairs(samples_array, pars=c("R","alpha")) #or with separate chains

# 95% central posterior intervals (CPI)
mean(samples$R)
mean(samples$alpha)
mean(samples$kD)
mean(samples$kE)
quantile(samples$R, c(0.025, 0.975))
quantile(samples$alpha, c(0.025, 0.975))
quantile(samples$kD, c(0.025, 0.975))
quantile(samples$kE, c(0.025, 0.975))

# R: Prior (weakly informative)
R_grid <- seq(0, 100, length.out=1000)
R_prior_dens <- dnorm(log(R_grid), mean=log(1), sd=sqrt(1/0.666))
plot(R_grid, R_prior_dens, type="l", ylim=c(0, max(R_prior_dens)))
plot(R_grid, R_prior_dens, type="l", ylim=c(0, max(R_prior_dens)), xlim=c(0,10))

# R: Posterior compared to prior
R_post_dens <- density(samples$R)
plot(R_post_dens, xlim=c(0, 100), col="blue")
lines(R_grid, R_prior_dens)
legend("topright", c("Prior","Posterior"), lty=1, col=c("black","blue"))

# alpha: Prior
alpha_grid <- seq(0, 10, length.out=100000)
alpha_prior_dens <- dnorm(log(alpha_grid), mean=log(0.001), sd=sqrt(1/0.5))
plot(alpha_grid, alpha_prior_dens, type="l", ylim=c(0, max(alpha_prior_dens)))
plot(alpha_grid, alpha_prior_dens, type="l", ylim=c(0, max(alpha_prior_dens)), xlim=c(0,0.1))
plot(alpha_grid, alpha_prior_dens, type="l", ylim=c(0, max(alpha_prior_dens)), xlim=c(0,0.01))

# alpha: Posterior compared to prior
alpha_post_dens <- density(samples$alpha)
plot(alpha_post_dens, xlim=c(0, 1), col="blue")
lines(alpha_grid, alpha_prior_dens)
legend("topright", c("Prior","Posterior"), lty=1, col=c("black","blue"))

# Posterior predictive checks
# Simulate the model using posterior samples for parameter values

# NBBG model
Ricker_NBBG <- function(Nt, R, alpha, kD, kE, p=0.5) {

    # Environmental stochasticity
    Rtx <- rgamma(length(Nt), shape=kE, scale=R/kE)

    # Heterogeneity in sex (binomial)
    females <- rbinom(length(Nt), Nt, p)

    # Heterogeneity in individual birth rate plus density-independent
    # survival
    births <- rnbinom( length(Nt), size = kD * females + (females==0),
            mu = (females / p) * Rtx )

    # Density dependent survival
    survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )

    return(survivors)
}

n_samples <- nrow(samples)
Ntp1_sim <- matrix(NA, nrow=n_samples, ncol=n)
for ( i in 1:n_samples ) {
    Ntp1_sim[i,] <- Ricker_NBBG(Nt, samples[i,"R"], samples[i,"alpha"], 
                                    samples[i,"kD"], samples[i,"kE"])
}

# Each row is a simulated dataset

# Simulation compared to real data is pretty good!
par(mfrow = c(3, 3))
sample_i <- sample(1:n_samples, 9)
x_max <- max(Ntp1_sim[sample_i,])
for ( i in sample_i ) {
    plot(Nt, Ntp1, pch = 16, xlab = "N_t", ylab = "N_{t+1}",
        ylim= c(0, x_max), main = paste("Sim", i))
    points(Nt, Ntp1_sim[i,], col = "red", pch = 1)
}
for ( i in sample_i ) {
    plot(Nt, Ntp1, pch = 16, xlab = "N_t", ylab = "N_{t+1}", xlim=c(0,150),
        main = paste("Sim", i))
    points(Nt, Ntp1_sim[i,], col = "red", pch = 1)
}
par(mfrow = c(1, 1)) #set it back

# Generic checks also look pretty good
ppc_dens_overlay(y = Ntp1, yrep = Ntp1_sim[1:200,])
ppc_stat(y = Ntp1, yrep = Ntp1_sim, stat = "mean")
ppc_stat(y = Ntp1, yrep = Ntp1_sim, stat = "sd")



