set.seed(8562)

# True parameters
imm <- 50
s <- 0.6
sigma_proc <- 0.1  #ln scale
sigma_obs  <- 0.1  #ln scale

T <- 50

N <- rep(NA, T)
D <- rep(NA, T)

# Initial abundance
N[1] <- 20

# Process model
for ( t in 1:(T-1) ) {
    mu <- imm + s * N[t]
    N[t+1] <- exp(rnorm(1, log(mu), sigma_proc))
}

# Observation model
for ( t in 1:T ) {
    D[t] <- exp(rnorm(1, log(N[t]), sigma_obs))
}

# Plot the simulated data
plot(N, type="l", lwd=2, ylim=range(c(0,N,D)))
points(D, pch=16, col="red")
legend("bottomright", c("True state","Observed"),
    lty=c(1,NA), pch=c(NA,16), col=c("black","red"))


# Bayesian training and inference with Stan

library(cmdstanr)
library(bayesplot)
library(posterior)

mod_obs  <- cmdstan_model("sinkpop_obs.stan")
mod_proc <- cmdstan_model("sinkpop_proc.stan")
mod_ss   <- cmdstan_model("sinkpop_ss.stan")

stan_data <- list(
    T = length(D),
    D = D
)

# Observation uncertainty model

init_fun <- function() {
  list(
    lnN1 = log(D[1]),
    imm = 50,
    s = 0.6,
    sigma_obs = 0.1
  )
}

set.seed(1257)
fit_obs <- mod_obs$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
)

fit_obs$summary(
    variables = c("lnN1", "imm", "s", "sigma_obs")
)

draws_N_obs <- as_draws_matrix(fit_obs$draws("N"))
N_hat_obs <- colMeans(draws_N_obs)     #estimated states

plot(1:T, D, pch=16, col="red", xlab="Time", ylab="Population size",
    ylim=c(0,200), main="Observation uncertainty model")
lines(1:T, N, lwd=2)
lines(1:T, N_hat_obs, col="blue", lty=2)
legend("topleft", c("True state","Observed","Estimated state"),
    lty=c(1, NA, 2), pch=c(NA, 16, NA), col=c("black", "red", "blue"))

# Process uncertainty model

init_fun <- function() {
  list(
    imm = 50,
    s = 0.6,
    sigma_proc = 0.1
  )
}

set.seed(292)
fit_proc <- mod_proc$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
)

fit_proc$summary(
    variables = c("imm", "s", "sigma_proc")
)

draws_lnNtp1_proc <- as_draws_matrix(fit_proc$draws("lnNtp1"))
N_hat_proc <- exp(colMeans(draws_lnNtp1_proc))     #estimated states

plot(1:T, D, pch=16, col="red", xlab="Time", ylab="Population size",
    ylim=c(0,200), main="Process uncertainty model")
lines(1:T, N, lwd=2)
points(2:T, N_hat_proc, col="blue", lty=2)
segments(1:(T-1), D[1:(T-1)], 2:T, N_hat_proc, col="blue", lty=2)
legend("topleft", c("True state","Observed","Estimated state"),
    lty=c(1, NA, 2), pch=c(NA, 16, NA), col=c("black", "red", "blue"))

# State space model (process + observation uncertainty)

init_fun <- function() {
  list(
    lnN1 = log(D[1]),
    imm = 50,
    s = 0.6,
    sigma_proc = 0.1,
    sigma_obs = 0.1,
    z = rep(0, T-1)
  )
}

set.seed(6455)
fit_ss <- mod_ss$sample(
    data = stan_data,
    init = init_fun,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.99,
    max_treedepth = 15
)

fit_ss$summary(
    variables = c("imm", "s", "sigma_proc", "sigma_obs")
)

summarize_draws(fit_ss$draws(), "mean", "sd", "rhat", "ess_bulk", "ess_tail") |>
    print(n=120)

samples_lnN_ss <- as_draws_matrix(fit_ss$draws("lnN"))
N_hat_ss <- exp(colMeans(samples_lnN_ss))     #estimated states

plot(1:T, D, pch=16, col="red", xlab="Time", ylab="Population size",
    ylim=c(0,200), main="Process + observation uncertainty model")
lines(1:T, N, lwd=2)
lines(1:T, N_hat_ss, col="blue", lty=2)
legend("topleft", c("True state","Observed","Estimated state"),
    lty=c(1, NA, 2), pch=c(NA, 16, NA), col=c("black", "red", "blue"))


samples_array <- fit_ss$draws(variables=c("imm", "s", "sigma_proc", "sigma_obs"))

mcmc_trace(samples_array, pars = "imm")
mcmc_trace(samples_array, pars = "s")
mcmc_trace(samples_array, "imm", window=c(0,500))
mcmc_trace(samples_array, "s", window=c(0,500))

mcmc_trace(samples_array, pars = "sigma_proc")
mcmc_trace(samples_array, pars = "sigma_obs")
mcmc_trace(samples_array, "sigma_proc", window=c(0,500))
mcmc_trace(samples_array, "sigma_obs", window=c(0,500))

# Posterior distributions
samples <- as_draws_df(samples_array)
mcmc_hist(samples, pars=c("imm","s","sigma_proc","sigma_obs"))
mcmc_pairs(samples, pars=c("imm","s","sigma_proc","sigma_obs"))

# Uncertainty intervals for state N

samples_lnN <- as_draws_df(fit_ss$draws("lnN"))
mcmc_hist(samples_lnN)
samples_N <- exp(as.matrix(samples_lnN))[,1:T] #trim extraneous columns
ci95 <- NULL
for ( k in colnames(samples_N) ) {
    ci_k <- quantile(samples_N[,k], c(0.025, 0.975))
    ci95 <- rbind(ci95, ci_k)
}

plot(1:T, D, pch=16, col="red", xlab="Time", ylab="Population size",
    ylim=c(0,200), main="Process + observation uncertainty model")
lines(1:T, N, lwd=2)
lines(1:T, N_hat_ss, col="blue", lty=2)
legend("topleft", c("True state","Observed","Estimated state"),
    lty=c(1, NA, 2), pch=c(NA, 16, NA), col=c("black", "red", "blue"))
# Add credible interval for N
lower <- ci95[,1]
upper <- ci95[,2]
polygon(
    x = c(1:T, T:1),
    y = c(lower, rev(upper)),
    col = adjustcolor("blue", alpha.f = 0.2),
    border = NA
)
