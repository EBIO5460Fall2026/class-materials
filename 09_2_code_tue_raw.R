library(R2jags)
library(bayesplot)

tribdata <- read.csv("data/ricker_data.csv")
Nt <- round(tribdata$At)
Ntp1 <- round(tribdata$Atp1)


# The JAGS/BUGS model definition
file.show("ricker_poisson.jags")

### MCMC: Rosenbluth algorithm (Metropolis-Hastings)

n = length(Nt)
data <- list(
    Nt = Nt,
    Ntp1 = Ntp1,
    n = n
)

params <- c("R", "alpha")

# Initial simulation to check chain behavior
Ricker_poisson_sim <- jags(
    data = data,
    parameters.to.save = params,
    model.file = "ricker_poisson.jags",
    n.chains = 1,
    n.iter = 4000,
    n.burnin = 0,      #keep all
    n.thin = 1         #iterations
)
Ricker_poisson_sim

samples <- as.data.frame(Ricker_poisson_sim$BUGSoutput$sims.matrix)
mcmc_trace(samples, "R")
mcmc_trace(samples, "R", window=c(0,500))
mcmc_trace(samples, "alpha")
mcmc_trace(samples, "alpha", window=c(0,500))
mcmc_acf(samples, c("R","alpha"))

# It doesn't look like much burn in is needed; 1000 should be conservative.
# There is no autocorrelation, so we don't need to aggressively thin.

# Multiple chains now and bumping up the iterations
Ricker_poisson_sim <- jags(
    data = data,
    parameters.to.save = params,
    model.file = "ricker_poisson.jags",
    n.chains = 4,
    n.iter = 6000,
    n.burnin = 1000,
    n.thin = 2
)
Ricker_poisson_sim

# Check convergence
samples_array <- Ricker_poisson_sim$BUGSoutput$sims.array
mcmc_trace(samples_array, "R")
mcmc_trace(samples_array, "R", window=c(0,500))
mcmc_trace(samples_array, "alpha", window=c(0,500))

# Posterior distributions (including 2D view)
samples <- as.data.frame(Ricker_poisson_sim$BUGSoutput$sims.matrix)
head(samples)
mcmc_hist(samples, pars=c("R","alpha"))
mcmc_pairs(samples, pars=c("R","alpha"))
mcmc_pairs(samples_array, pars=c("R","alpha")) #or with separate chains
with(samples, plot(alpha, R))

# 95% central posterior intervals (CPI)
quantile(samples$R, c(0.025, 0.975))
quantile(samples$alpha, c(0.025, 0.975))

# R: Prior (weakly informative)
precision <- 0.25
R_grid <- seq(0, 100, length.out=1000)
R_prior_dens <- dnorm(log(R_grid), mean=log(1), sd=sqrt(1/precision))
plot(R_grid, R_prior_dens, type="l", ylim=c(0, max(R_prior_dens)))
plot(R_grid, R_prior_dens, type="l", ylim=c(0, max(R_prior_dens)), xlim=c(0,10))

# R: Posterior compared to prior
R_post_dens <- density(samples$R)
plot(R_post_dens, xlim=c(0, 10), col="blue")
lines(R_grid, R_prior_dens)
legend("topright", c("Prior","Posterior"), lty=1, col=c("black","blue"))

# alpha: Prior
precision <- 0.25
alpha_grid <- seq(0, 10, length.out=100000)
alpha_prior_dens <- dnorm(log(alpha_grid), mean=log(0.001), sd=sqrt(1/precision))
plot(alpha_grid, alpha_prior_dens, type="l", ylim=c(0, max(alpha_prior_dens)))
plot(alpha_grid, alpha_prior_dens, type="l", ylim=c(0, max(alpha_prior_dens)), xlim=c(0,0.1))
plot(alpha_grid, alpha_prior_dens, type="l", ylim=c(0, max(alpha_prior_dens)), xlim=c(0,0.01))

# alpha: Posterior compared to prior
alpha_post_dens <- density(samples$alpha)
plot(alpha_post_dens, xlim=c(0, 0.01), col="blue")
lines(alpha_grid, alpha_prior_dens)
legend("topright", c("Prior","Posterior"), lty=1, col=c("black","blue"))

# Posterior predictive checks
# Simulate the model using posterior samples for parameter values

n_samples <- nrow(samples)
Ntp1_sim <- matrix(NA, nrow=n_samples, ncol=n)

Ricker_Poisson <- function(Nt, R, alpha) {

    # Births and density independent survival
    births_di_survive <- rpois( length(Nt), Nt * R )

    # Density dependent survival
    survivors <- rbinom( length(Nt), births_di_survive, exp( -1 * alpha * Nt ) )

    return(survivors)
}

for ( i in 1:n_samples ) {
    Ntp1_sim[i,] <- Ricker_Poisson(Nt, samples[i,"R"], samples[i,"alpha"])
}

# Each row is a simulated dataset

# The most telling plot is often direct data simulation compared to real data
par(mfrow = c(3, 3))
for ( i in sample(1:n_samples, 9) ) {
    plot(Nt, Ntp1, pch = 16, xlab = "N_t", ylab = "N_{t+1}",
        main = paste("Sim", i))
    points(Nt, Ntp1_sim[i,], col = "red", pch = 1)
}
#
for ( i in sample(1:n_samples, 9) ) {
    plot(Nt, Ntp1, pch = 16, xlab = "N_t", ylab = "N_{t+1}", xlim=c(0,150),
        main = paste("Sim", i))
    points(Nt, Ntp1_sim[i,], col = "red", pch = 1)
}
par(mfrow = c(1, 1)) #set it back

# Some generic checks also pick up issues
ppc_dens_overlay(y = Ntp1, yrep = Ntp1_sim[1:200,])
ppc_stat(y = Ntp1, yrep = Ntp1_sim, stat = "mean")
ppc_stat(y = Ntp1, yrep = Ntp1_sim, stat = "sd")
