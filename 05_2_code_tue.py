
import numpy as np

rng = np.random.default_rng(123)

b = 0.2
d = 0.1

alive = np.full(10, 1)
b = np.full(N_0, 0.25)    # birth rate
m = np.full(N_0, 0.2)     # death rate

t = 0
t_max = 30
N = sum(alive)
while t < t_max and N > 0:



    total_birth = np.sum(b * alive)
    total_death = np.sum(m * alive)
    intensity = total_birth + total_death

    # time until next event
    t += rng.exponential(scale=1/total_rate)

    # determine event type
    if rng.binomial(1, total_birth / intensity):

        # birth
        alive = np.append(alive, 1)

    else:

        # death
        i = rng.choice(np.where(alive == 1)[0])
        alive[i] = 0

    N = sum(alive)

    print(t, N)
