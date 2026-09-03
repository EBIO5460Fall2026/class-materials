
Useful distributions to get started with stochastic modeling

| Distribution | Type | R | SciPy | NumPy |
|--------------|--|----|--------|--------|
| Uniform | Continuous | unif | uniform | uniform |
| Bernoulli | Discrete | binom(size=1) | bernoulli | binomial(n=1) |
| Binomial | Discrete | binom | binom | binomial |
| Exponential | Continuous | exp | expon | exponential |
| Poisson | Discrete | pois | poisson | poisson |
| Negative Binomial | Discrete | nbinom | nbinom | negative_binomial |
| Normal | Continuous | norm | norm | normal |
| Gamma | Continuous | gamma | gamma | gamma |
| Beta | Continuous | beta | beta | beta |
| Lognormal | Continuous | lnorm | lognorm | lognormal |

Usage: Python
```python
import numpy as np
from scipy import stats

rng = np.random.default_rng()
np.rng.poisson()              # Random numbers
x = ...
stats.poisson.pmf(x, mu=5)    # PMF (discrete) or use .pdf for PDF (continuous)
stats.poisson.cdf(x, mu=5)    # CDF
stats.poisson.ppf(0.95, mu=5) # Quantiles
```

Usage: R
```R
set.seed(12345)
x = 0:25
rpois(100, lambda=5)     # Random numbers
dpois(x, lambda=5)       # PMF/PDF
ppois(x, lambda=5)       # CDF
qpois(0.95, lambda=5)    # Quantiles
```