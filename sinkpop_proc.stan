data {
    int<lower=1> T;
    vector[T] D;
}

parameters {
    real<lower=0> imm;
    real<lower=0,upper=1> s;
    real<lower=0> sigma_proc;
}

transformed parameters {
    vector[T-1] lnNtp1;
    for (t in 1:(T-1)) {
        real mu = imm + s * D[t];
        lnNtp1[t] = log(mu);
    }
}

model {
    // Priors
    imm ~ normal(0, 100);
    s ~ beta(2, 2); //middle values
    sigma_proc ~ normal(0, 0.1);
    
    // Likelihood
    D[2:T] ~ lognormal(lnNtp1, sigma_proc);
}
