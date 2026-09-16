# Continuous-time birth-death IBM

import numpy as np
import matplotlib.pyplot as plt

# Simulation control
rng = np.random.default_rng(95756)
t = 0
t_max = 30

# Initial biological state

# --individual scale
N_0 = 10
alive = np.full(N_0, 1)
b = np.full(N_0, 0.25)    # birth rate
m = np.full(N_0, 0.2)     # death rate

# --population scale
N_t = np.sum(alive) #scale transition (sum over individuals -> population)
N = N_t      #track N
times = t    #track time

# Simulate events
while t < t_max and N_t > 0:

    # calculate event intensity
    total_birth = np.sum(b * alive)
    total_death = np.sum(m * alive)
    intensity = total_birth + total_death

    # waiting time to next event
    w = rng.exponential(scale=1/intensity)

    # advance time
    t = t + w

    # which event type occurs?
    event_birth = rng.binomial(1, total_birth / intensity) #1 = birth

    # which individual does the event happen to?
    i = rng.choice(np.where(alive == 1)[0])  #equal p, identical individuals

    # update system state (do event)
    if event_birth:
        # add a new individual and its traits
        alive = np.append(alive, 1)
        b = np.append(b, 0.25)
        m = np.append(m, 0.2)
    else:
        # death
        alive[i] = 0

    N_t = sum(alive)
    N = np.append(N, N_t)
    times = np.append(times, t)

print(alive)

plt.figure()
plt.step(times, N, where='post')
plt.ylim(bottom=0)
plt.xlabel("Time")
plt.ylabel("N")
plt.title("Birth-death IBM")
