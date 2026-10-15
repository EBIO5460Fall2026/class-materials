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

    lnR     ~ normal(log(1), 2);
    lnalpha ~ normal(log(0.001), 2);
    lnkD    ~ normal(log(1), 2);
    lnkE    ~ normal(log(1), 2);

    for (i in 1:n) {

        Rtx[i] ~ gamma(kE, kE / R);

        # log probability for a single observation
        vector[Nt[i] + 1] lp;

        real s = exp(-alpha * Nt[i]);

        for (f in 0:Nt[i]) {

            if (f == 0) {

                # vector indexing starts at 1, so index = f + 1
                lp[f+1] = binomial_lpmf(0 | Nt[i], 0.5);

                if (Ntp1[i] == 0)
                    # p = 1; (ln(p) = 0)
                    lp[f+1] += 0;
                else
                    # p = 0; (ln(p) = -Inf)
                    lp[f+1] += negative_infinity();

            } else {

                real r = kD * f;
                real mu = (f / 0.5) * Rtx[i] * s;
                lp[f+1] =
                    binomial_lpmf(f | Nt[i], 0.5) +
                    neg_binomial_2_lpmf(Ntp1[i] | r, mu);
            }
        }

        target += log_sum_exp(lp);  #log posterior (Stan keyword)
    }
}
