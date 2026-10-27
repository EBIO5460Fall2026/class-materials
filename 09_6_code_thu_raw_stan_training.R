
# Setup:
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
# library(cmdstanr)
# check_cmdstan_toolchain(fix = TRUE)
# install_cmdstan()

# CmdStanR will work in parallel by default. To see how many cores you have:
parallel::detectCores(logical = FALSE)

# Tools
library(cmdstanr)
library(bayesplot)
library(posterior)
library(coda) #for HPDI

# Read data
tribdata <- read.csv("data/ricker_data.csv")
Nt <- round(tribdata$At)
Ntp1 <- round(tribdata$Atp1)

#----Ricker PB model-------------------------------
# We'll start with a simpler model of demographic stochasticity + sex

n <- length(Nt)
data <- list(
    Nt = Nt,
    Ntp1 = Ntp1,
    n = n
)

# Compile model
mod <- cmdstan_model("ricker_pb.stan")

# Useful here to monitor CPUs on your computer to confirm 4 cores in use

Ricker_pb_sim <- mod$sample(
    data = data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 4000,
    seed = 5253
)

samples_array <- Ricker_pb_sim$draws(variables = c("R", "alpha"))
summarise_draws(samples_array, mean, sd, rhat, ess_bulk, ess_tail, mcse_mean)

mcmc_trace(samples_array, pars = "R")
mcmc_trace(samples_array, pars = "alpha")
mcmc_trace(samples_array, "R", window=c(0,500))
mcmc_trace(samples_array, "alpha", window=c(0,500))

# Posterior distributions
samples <- as_draws_df(samples_array)
head(samples)
mcmc_hist(samples, pars=c("R","alpha"))
mcmc_pairs(samples, pars=c("R","alpha"))

# 95% central posterior intervals (CPI)
quantile(samples$R, c(0.025, 0.975))
quantile(samples$alpha, c(0.025, 0.975))


#---NBBG model---------------------------------
# Now the most complex model

n <- length(Nt)
data <- list(
    n = n,
    Nt = Nt,
    Ntp1 = Ntp1
)

# Compile model
mod <- cmdstan_model("ricker_nbbg.stan")

# Useful here to monitor CPUs on your computer to confirm 4 cores in use

Ricker_nbbg_sim <- mod$sample(
    data = data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 4000,
    seed = 5253
)

samples_array <- Ricker_nbbg_sim$draws(c("R", "alpha", "kD", "kE"))
summarise_draws(samples_array, mean, sd, rhat, ess_bulk, ess_tail, mcse_mean)

mcmc_trace(samples_array, pars="R")
mcmc_trace(samples_array, pars="alpha")
mcmc_trace(samples_array, "R", window=c(0,500))
mcmc_trace(samples_array, "alpha", window=c(0,500))
mcmc_trace(samples_array, pars="kD")
mcmc_trace(samples_array, pars="kE")
mcmc_trace(samples_array, "kD", window=c(0,500))
mcmc_trace(samples_array, "kE", window=c(0,500))
mcmc_trace(samples_array, pars="kD", transformations="log")
mcmc_trace(samples_array, pars="kE", transformations="log")
mcmc_acf(samples_array, c("R","alpha","kD","kE"))

# Posterior distributions
samples <- as_draws_df(samples_array)
mcmc_hist(samples, pars=c("R","alpha","kD","kE"))
mcmc_hist(samples, pars="kD") +
    xlim(0, 50)
mcmc_pairs(samples, pars=c("R","alpha","kD","kE"))

# 95% central posterior intervals (CPI)
quantile(samples$R, c(0.025, 0.975))
quantile(samples$alpha, c(0.025, 0.975))
quantile(samples$kD, c(0.025, 0.975))
quantile(samples$kE, c(0.025, 0.975))

# For the k parameters it makes sense to use highest posterior density intervals
# (HPDI) since they are very asymmetric
HPDinterval(as.mcmc(samples$kD))
HPDinterval(as.mcmc(samples$kE))


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

samples <- as.data.frame(samples)
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



# LOOIC
# LOOCV where predictive performance is measured by the posterior predictive
# probability of the test observation (held out data point), averaged over the
# posterior samples from the model trained on the remaining observations.

# Helper functions to aggregate exponents of log probabilities that might be
# large negative numbers (as log probabilities often are). Prevents underflow of
# exp(lp), i.e. prevents exp(lp) = 0 when lp is << 0.

log_sum_exp <- function(x) {
    m <- max(x)
    return(m + log(sum(exp(x - m))))
}

log_mean_exp <- function(x) {
    m <- max(x)
    return(m + log(mean(exp(x - m))))
}

# Generate one sample of the log predictive probability for a test observation
# derived from a sample of the posterior. In the Bayesian model-comparison
# literature this quantity is commonly referred to as a log predictive density
# (lpd), even for discrete data. This is a vectorized R version of the same code
# in the model block of the Stan code.
#
# Ntp1, Nt:         the held-out (test) observation
# R, alpha, kD, kE: one sample from the posterior of the model trained without
#                   the test observation.
#  
sample_lpd <- function(Ntp1, Nt, R, alpha, kD, kE) {
    p_fem <- 0.5
    Rtx <- rgamma(1, kE, kE / R) #random draw
    s <- exp(-alpha * Nt)
    f <- 0:Nt
    mu <- (Rtx / p_fem) * f * s
    phi <- kD * f + (f==0) #adds 1 when females=0 to avoid NaN's
    lpd <- dbinom(f, Nt, p_fem, log=TRUE) + 
           dnbinom(Ntp1, size=phi, mu=mu, log=TRUE)
    return(log_sum_exp(lpd))    
}

# Train the model with observation i left out and return the log posterior
# predictive probability of the held-out observation. Later, summing these
# values across held-out observations gives the ELPD estimate used to compute
# LOOIC.
#
loo_lpd_i <- function(i, Ntp1, Nt) {

    datalist = list(
        n = length(Ntp1) - 1,
        Nt = Nt[-i],
        Ntp1 = Ntp1[-i]
    )

    fit <- mod$sample(
        data = datalist,
        chains = 1,
        iter_warmup = 1000,
        iter_sampling = 1000
    )

    samples <- as_draws_df(fit$draws(c("R", "alpha", "kD", "kE")))
    n_samples <- nrow(samples)
    
    lpd_samples <- rep(NA, n_samples)
    for ( j in 1:n_samples ) {
        lpd_samples[j] <-
            sample_lpd(
                Ntp1[i],
                Nt[i],
                samples$R[j],
                samples$alpha[j],
                samples$kD[j],
                samples$kE[j]
            )
    }

    # Monte Carlo estimate of the log posterior predictive probability
    return(log_mean_exp(lpd_samples))
}


# Parallel processing
# We can speed this up using parallel computing. To use parallel `for` loops,
# we'll use the `doFuture` package.
# 
library(doFuture) #For parallel processing
library(doRNG) #For safe random numbers in parallel processing
registerDoFuture()


# To prepare for a parallel run, we first set a "plan" for the parallel
# workflow, which I'm here indicating should have 12 workers, using 12 of the 16
# processor cores I have available on my laptop. This will set up 12 separate R
# environments to run in parallel. You can omit the `workers` argument, which
# will default to the number of available cores on your computer.

plan(multisession, workers=12)

# Next we do: "for each observation i, do loo_lpd_i". The `%dorng%` operator
# ensures that random number generation is done properly and is reproducible in
# parallel (this is done by the `doRNG` package, loaded above). Combine says to
# make a vector (c), and we need to pass any packages needed by the loo_lpd_i
# function to the workers. 

set.seed(90076)
lpd_i <- foreach(i = 1:n,
    .combine = c,
    .packages = c("cmdstanr", "posterior")
) %dorng% {
    loo_lpd_i(i, Ntp1, Nt)
}

# save(lpd_i, file="09_6_lpdi_saved.Rdata")
load("09_6_lpdi_saved.Rdata")

elpd_loo <- sum(lpd_i)
looic <- -2 * elpd_loo

se_elpd <- sqrt(length(lpd_i) * var(lpd_i))
se_looic <- 2 * se_elpd

# We can now use this to compare models (see refresher topics)
print(paste("LOOIC: ", looic, "+/- se:", se_looic))
