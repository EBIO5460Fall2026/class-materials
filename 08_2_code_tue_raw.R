
Ricker_Poisson_base <- function(Nt, beta, m, alpha){
    # Births
    births <- rpois( length(Nt), Nt * beta )
    # Density dependent and density independent survival
    survivors <- rbinom( length(Nt), births, (1 - m) * exp( -1 * alpha * Nt ) )
    return(survivors)
}

# It's not possible to distinguish between beta and m from data. This is called
# non-identifiability. We can change one or the other for identical effect, not
# just in the expectation (mean of the stochastic model) but also in the full
# distribution. This is a consequence of how the different stochastic processes
# combine together in this case. They form a composite parameter, which is the
# finite growth rate of the population:
#
#     R = beta ( 1 - m )
#
# Poisson births, Pois(beta), combined with binomial, density-independent
# survival, Binom(births, 1-m), gives Poisson survivors.

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

# Then we account for density-dependent survival. This is less biologically
# intuitive because density-independent survival and density-dependent survival
# actually happen at the same time, not one after the other. The model, while
# correct in outcome, is not following the biological algorithm. But this allows
# us to combine the non-identifiable parameters.

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


# Combining this into a function

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

Ricker_NBBG <- function(Nt, R, alpha, kD, kE, p=0.5) {

    # Heterogeneity in R (birth rate/DI mortality) between times or locations
    # (gamma with mean R)
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



# Vectorized simulation
sims <- 1000000
Nt <- rep(25, sims)
R <- 5 * ( 1 - 0.2 )
alpha <- 0.002
kD <- 10
kE <- 10
p <- 0.5
Ntp1 <- Ricker_NBBG(Nt, R, alpha, kD, kE, p)

# Distribution (PMF) of Ntp1
hist(Ntp1, breaks=seq(min(Ntp1)-0.5, max(Ntp1)+0.5, by=1),
    probability=TRUE, col="lightblue", xlab="number surviving",
    ylab="Probability mass", main="")


# Show that the theoretical PMF is correct
source("source/Ricker_dnbinombinomgamma.R")

# Add theoretical PMF to plot
X <- min(Ntp1):max(Ntp1)
pmf <- Ricker_dnbinombinomgamma(X, rep(Nt, length(X)), p, R, alpha, kD, kE)
points(X, pmf, col="darkorange")
title("NBBG(Nt R exp(-alpha Nt)) (points) vs simulation (bars)")


#----Likelihood for the Poisson Ricker model
#
Ricker_pois.nll <- function(lnR, lnalpha){
  #Probability of the data
  mu <- Nt * exp(lnR) * exp(-exp(lnalpha) * Nt)
  -sum(dpois(Ntp1, mu, log=TRUE))
}

#----Likelihood for the Negative binomial-binomial-gamma Ricker model
#
Ricker_nbinombinomgamma.nll <- function(lnR, lnalpha, lnkD, lnkE){
  #Probability of the data
# source("source/Ricker_var.R
  -sum(Ricker_dnbinombinomgamma(Ntp1, Nt, 0.5, exp(lnR), exp(lnalpha),
                                          exp(lnkD), exp(lnkE), log=TRUE))
}


# source("source/Ricker_var.R") #Variances of stochastic Ricker models
library(stats4) #mle

tribdata <- read.csv("data/ricker_data.csv")
tribdata$Nt <- round(tribdata$At)
tribdata$Ntp1 <- round(tribdata$Atp1)
attach(tribdata)


# Fit Poisson Ricker

llfit <- mle( Ricker_pois.nll,
              start=list(
                        lnR=log(2.613),
                        lnalpha=log(0.003731)
                    )
            )
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha")
phat



# Fit Negative binomial-binomial-gamma Ricker

llfit <- mle( Ricker_nbinombinomgamma.nll,start=list(
                                             lnR=log(2.613),
                                             lnalpha=log(0.003731),
                                             lnkD=log(1.1475),
                                             lnkE=log(26.6221) ) )
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat



# Figure

# Need the basic Ricker function here

m_Nt <- expression( italic(N) [ italic(t) ] )
m_Ntp1 <- expression( italic(N) [ italic(t)+1 ] )
plot(Nt,Ntp1,xlab=m_Nt,ylab=m_Ntp1,col="black")
lines(1:1100,Ricker(1:1100,phat["R"], phat["alpha"]),col="red") #col="grey75"
x <- c(2,6,10,18,28,46,78,130,215,360,600,1000)
y <- Ricker(x,phat["R"], phat["alpha"])
# std <- sqrt(Ricker_nbinombinomgamma.var(x,phat["R"], phat["alpha"],
#                                           phat["kD"],phat["kE"]))
segments(x,y-std,x,y+std,col="red") #col="grey75"
