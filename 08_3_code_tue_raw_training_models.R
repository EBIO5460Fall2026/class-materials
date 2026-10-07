# Confronting models with data
# Maximum likelihood training algorithm
# We are estimating the model parameters from experimental data

source("source/Ricker_dnbinombinomgamma.R")
library(stats4) #mle

tribdata <- read.csv("data/ricker_data.csv")
tribdata$Nt <- round(tribdata$At)
tribdata$Ntp1 <- round(tribdata$Atp1)
attach(tribdata)

# First we need likelihood functions
# The set up here embeds the data, N_{t} and N_{t+1}, as global variables in the
# likelihood function. The parameters are all continous positive quantities, so
# we pass them to the likelihood functions on the log scale. This allows the
# optimizer algorithm to range freely from -Inf to +Inf. The likelihood
# functions return the negative log likelihood of the data.

# Likelihood for the Poisson Ricker model
Ricker_pois.nll <- function(lnR, lnalpha) {
    #Probability of the data
    mu <- Nt * exp(lnR) * exp(-exp(lnalpha) * Nt)
    -sum(dpois(Ntp1, mu, log=TRUE))
}

# Likelihood for the Negative binomial-binomial-gamma Ricker model
Ricker_nbinombinomgamma.nll <- function(lnR, lnalpha, lnkD, lnkE) {
    #Probability of the data
    -sum(Ricker_dnbinombinomgamma(Ntp1, Nt, 0.5, exp(lnR), exp(lnalpha),
                                          exp(lnkD), exp(lnkE), log=TRUE))
}

# Now we can optimize the likelihoods - find the values of the parameters that
# minimize the negative log likelihood, hence maximize the likelihood. Here we
# use the function mle for maximum likelihood estimation. This is essentially a
# wrapper function to the more fundamental optim() function in base R. The
# default training algorithm is a robust descent algorithm called "BFGS". We
# need pretty good starting values for the parameters, and it's a good idea to
# try different starting values to ensure we find the optimum. After training,
# we compare the models with AIC, an approximation to the leave one out cross
# validation algorithm (LOOCV).

# Fit Poisson Ricker

llfit <- mle( Ricker_pois.nll,
              start=list(
                        lnR=log(2.5),
                        lnalpha=log(0.004)
                    )
            )
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha")
phat
AIC(llfit)

# Fit Negative binomial-binomial-gamma Ricker

llfit <- mle( Ricker_nbinombinomgamma.nll,start=list(
                                             lnR=log(2.5),
                                             lnalpha=log(0.004),
                                             lnkD=log(1),
                                             lnkE=log(25) ) )
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat
AIC(llfit)


# Plot the best fitting model (NBBG Ricker) with the data

# The mean of the NBBG model is the deterministic Ricker model
Ricker_mu <- function(Nt, R, alpha) {
    mu = Nt * R * exp(-1 * alpha * Nt)
    return(mu)
}

m_Nt <- expression( italic(N) [ italic(t) ] )
m_Ntp1 <- expression( italic(N) [ italic(t)+1 ] )
plot(Nt, Ntp1, xlab=m_Nt, ylab=m_Ntp1, col="black")
lines(1:1100, Ricker_mu(1:1100, phat["R"], phat["alpha"]), col="red")

# As a simple representation of the stochastic component, add the standard
# deviation of the NBBG model. A function for the variance is with the PMF in
# the source directory.
x <- c(2,6,10,18,28,46,78,130,215,360,600,1000)
y <- Ricker_mu(x, phat["R"], phat["alpha"])
std <- sqrt(Ricker_nbinombinomgamma.var(x, phat["R"], phat["alpha"],
                                           phat["kD"],phat["kE"]))
segments(x, y - std, x, y + std, col="red")
