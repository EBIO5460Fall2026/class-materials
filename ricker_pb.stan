data {
    int<lower=1> n;
    array[n] int<lower=0> Nt;
    array[n] int<lower=0> Ntp1;
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

    lnR     ~ normal(log(1), 1.2);
    lnalpha ~ normal(log(0.001), 1.4);

    real p_fem = 0.5;

    for (i in 1:n) {

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

                real mu = (R / p_fem) * f * s;
                lp_f[f+1] =
                    binomial_lpmf(f | Nt[i], p_fem) +
                    poisson_lpmf(Ntp1[i] | mu);
            }

        }

        target += log_sum_exp(lp_f);  //log posterior (Stan keyword)
    }
}

// Version that doesn't work
// It probably should but it consistently has problems with initial values

// model {

//     lnR     ~ normal(log(1), 2);
//     lnalpha ~ normal(log(0.001), 2);

//     for (i in 1:n) {

//         // log probability for each possible value of f
//         vector[Nt[i] + 1] lp_f;

//         real s = exp(-alpha * Nt[i]);

//         for (f in 0:Nt[i]) {
//             real mu = (R / 0.5) * f * s;
//             lp_f[f+1] =
//                 binomial_lpmf(f | Nt[i], 0.5) +
//                 poisson_lpmf(Ntp1[i] | mu);
//         }

//         target += log_sum_exp(lp_f);  //log posterior (Stan keyword)
//     }
// }

