data {
    int<lower=1> T;
    array[T] int<lower=0> Nt;
    array[T] int<lower=0> Ntp1;
}

parameters {
    real lnR;
    real lnalpha;
}

transformed parameters {
    real<lower=0> R = exp(lnR);
    real<lower=0> alpha = exp(lnalpha);
}

model {

    // example priors
    lnR ~ normal(0, 2);
    lnalpha ~ normal(-5, 2);

    for(t in 1:T){
        real mu;
        mu = Nt[t] * R * exp(-alpha * Nt[t]);
       Ntp1[t] ~ poisson(mu);
    }
}
