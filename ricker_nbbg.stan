data {
    int<lower=1> n;
    array[n] int<lower=0> Nt;
    array[n] int<lower=0> Ntp1;
}

parameters {
    real lnR;
    real lnalpha;
    real lnkD;
    real lnkE;

    vector<lower=0>[n] Rtx;
}

transformed parameters {
    real<lower=0> R = exp(lnR);
    real<lower=0> alpha = exp(lnalpha);
    real<lower=0> kD = exp(lnkD);
    real<lower=0> kE = exp(lnkE);
}

model {

    lnR     ~ normal(log(1), 1.2);
    lnalpha ~ normal(log(0.001), 1.5);
    lnkD    ~ normal(log(10), 1.5);
    lnkE    ~ normal(log(10), 1.5);

    real p_fem = 0.5;

    for (i in 1:n) {

        Rtx[i] ~ gamma(kE, kE / R);

        // log probability for each possible value of f
        vector[Nt[i] + 1] lp_f;

        real s = exp(-alpha * Nt[i]);

        for (f in 0:Nt[i]) {

            if (f == 0) {

                // handle f = 0 explicitly
                // vector indexing starts at 1, so index = f + 1
                if (Ntp1[i] == 0)
                    // p = dbinom(0, 0.5)
                    lp_f[f+1] = binomial_lpmf(0 | Nt[i], p_fem);
                else
                    // p = 0; (ln(p) = -Inf)
                    lp_f[f+1] = negative_infinity();

            } else {

                real phi = kD * f;    //dispersion depends on f
                real mu = (Rtx[i] / p_fem) * f * s;
                lp_f[f+1] =
                    binomial_lpmf(f | Nt[i], p_fem) +
                    neg_binomial_2_lpmf(Ntp1[i] | mu, phi);
            }
        }

        // The total probability across f values is the sum of the probabilities
        // (i.e events f = 0 OR f = 1 OR f = 2 ...). This function sums the
        // probabilities (p_f = exp(lp_f)) then takes the log of the sum.
        target += log_sum_exp(lp_f);  //log posterior (Stan keyword)
    }
}
