# Exponential distribution and event waiting times

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom

# Initiate the random number generator

rng = np.random.default_rng(2293)

c = 0.02         # roughly probability of death per year (so, lives ca 50 years)
delta_t = 1/365  # time increment 1 day, units are years
p = c * delta_t  # probability per day
print(p)

state = 0
t = 0
while state == 0:
    t = t + delta_t
    if ( rng.binomial(1, p, 1) ):
        state = 1
print(t)


# Turn this algorithm into a function

def death_continous_time(c, delta_t):
    p = c * delta_t
    alive = 1
    t = 0
    while alive:
        t = t + delta_t
        if ( rng.binomial(1, p, 1) ):
            alive = 0
    return t


# Conduct a simulation study

rng = np.random.default_rng(58773)

c = 0.02
delta_t = 1/365
reps = 1000
time_to_death = np.empty(reps)

for i in range(reps):
    time_to_death[i] = death_continous_time(c, delta_t)
    if i % 100 == 0:
        print(i)

print(time_to_death[:20])

plt.figure()
plt.hist(time_to_death, bins=20, density=True)
plt.xlabel("Time to death")
plt.ylabel("Density")
plt.show()

# Theoretical exponential distribution

t = np.linspace(0, np.max(time_to_death) + 1, 100)
pdf = c * np.exp(-c * t)
plt.plot(t, pdf)

print(np.mean(time_to_death)) #Expected: 1 / c
print(1 / c**2) #Expected variance: 1 / c^2
print(np.var(time_to_death, ddof=1))
print(np.std(time_to_death, ddof=1)) # Expected standard deviation = mean



# Poisson process simulation
#
# Simulate a sequence of death events through time.
#
# The waiting time until the next death is exponentially distributed.
# After each death occurs, we draw another waiting time and repeat.

rng = np.random.default_rng(481)

rate = 0.5          # deaths per unit time
n_events = 100

# Store results

waiting_time = np.empty(n_events)
death_time = np.empty(n_events)
deaths = np.arange(1, n_events + 1)

# First death

waiting_time[0] = rng.exponential(scale=1/rate)
death_time[0] = waiting_time[0]

# Remaining deaths

for i in range(1, n_events):
    waiting_time[i] = rng.exponential(scale=1/rate)
    death_time[i] = death_time[i-1] + waiting_time[i]

print(waiting_time)
print(death_time)



plt.figure()

plt.step(death_time, deaths, where='post')

plt.xlabel("Time")
plt.ylabel("Cumulative deaths")

plt.title("Poisson process")
plt.show()


