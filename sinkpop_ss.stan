data {
    int<lower=1> T;
    vector<lower=0>[T] D;
}

parameters {
    real lnN1;
    vector[T-1] z;
    real<lower=0> imm;
    real<lower=0,upper=1> s;
    real<lower=0> sigma_proc;
    real<lower=0> sigma_obs;
}

transformed parameters {
    vector[T] lnN;
    lnN[1] = lnN1;
    for(t in 1:(T-1)) {
        real mu = imm + s * exp(lnN[t]);
        lnN[t+1] = log(mu) + sigma_proc * z[t];
    }
}

model {
    // Priors
    lnN1 ~ normal(log(D[1]), 0.3);
    imm ~ normal(0, 100);
    s ~ beta(2, 2); //middle values
    // s ~ beta(8, 2); //mean 0.8
    //s ~ beta(1, 1); //no preference (uniform?)
    sigma_proc ~ normal(0, 0.1);
    sigma_obs ~ normal(0, 0.1);

    // Likelihood
    z ~ normal(0, 1);
    D ~ lognormal(lnN, sigma_obs);
}
