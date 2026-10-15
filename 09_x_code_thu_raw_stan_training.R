
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

# Read data
tribdata <- read.csv("data/ricker_data.csv")
Nt <- round(tribdata$At)
Ntp1 <- round(tribdata$Atp1)


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

mcmc_trace(samples_array, pars = "R")
mcmc_trace(samples_array, pars = "alpha")

summarise_draws(samples_array, mean, sd, rhat, ess_bulk, ess_tail, mcse_mean)


samples <- as_draws_df(samples_array)


