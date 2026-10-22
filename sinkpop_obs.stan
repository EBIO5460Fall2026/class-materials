data {
    int<lower=1> T;
    vector[T] D;
}

parameters {
    real<lower=0> lnN1;
    real<lower=0> imm;
    real<lower=0,upper=1> s;
    real<lower=0> sigma_obs;
}

transformed parameters {
    vector[T] N;
    N[1] = exp(lnN1);
    for(t in 1:(T-1)) {
        N[t+1] = imm + s * N[t];
    }
}

model {
    // Priors
    lnN1 ~ normal(log(D[1]), 0.3);
    imm ~ normal(0, 100);
    s ~ beta(2, 2); //middle values
    sigma_obs ~ normal(0, 0.1);
    
    // Likelihood
    D ~ lognormal(log(N), sigma_obs);
}
