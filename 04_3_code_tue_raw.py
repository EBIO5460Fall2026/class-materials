# Simulation: continuous time algorithm

# Exponential distribution and event waiting times

import numpy as np
import matplotlib.pyplot as plt

# Initiate the random number generator

rng = np.random.default_rng(2293)

m = 0.02         # roughly probability of death per year (so, lives ca 50 years)
delta_t = 1/365  # time increment 1 day, units are years
p = m * delta_t  # probability per day
print(p)

state = 0
t = 0
while state == 0:
    t = t + delta_t
    if ( rng.binomial(1, p, 1) ):
        state = 1
print(t)


# Turn this algorithm into a function

def death_continous_time(m, delta_t, rng):
    p = m * delta_t
    alive = 1
    t = 0
    while alive:
        t = t + delta_t
        if ( rng.binomial(1, p, 1) ):
            alive = 0
    return t


# Conduct a simulation study

rng2 = np.random.default_rng(58773)

m = 0.02
delta_t = 1/365
sims = 1000
time_to_death = np.empty(sims)

for i in range(sims):
    time_to_death[i] = death_continous_time(m, delta_t, rng2)
    if i % 100 == 0:
        print(i)

print(time_to_death[:20])

plt.figure()
plt.hist(time_to_death, bins=20, density=True)
plt.xlabel("Time to death")
plt.ylabel("Density")

# Theoretical exponential distribution

t = np.linspace(0, np.max(time_to_death) + 1, 100)
pdf = m * np.exp(-m * t)
plt.plot(t, pdf)

print(np.mean(time_to_death)) #Expected: 1 / m
print(1 / m**2) #Expected variance: 1 / m^2
print(np.var(time_to_death, ddof=1))
print(np.std(time_to_death, ddof=1)) # Expected standard deviation = mean

# Saving/reloading results from a long simulation

np.save("time_to_death.npy", time_to_death)
time_to_death = np.load("time_to_death.npy")

time_to_death_rl = np.load("time_to_death.npy")

#Same
time_to_death - time_to_death_rl



# Poisson process simulation
#
# Individual searching for prey
# Simulate prey discovery events up to time t_max.
#

rng = np.random.default_rng(1484)

s = 0.05          # search rate (ca probability per unit time)
t_max = 500       # stop simulation at this time
print(1/s)        # average time to find next prey item

# Store results
waiting_time = []
dinner_time = []

t = 0

while t < t_max:

    # draw waiting time until next prey found
    w = rng.exponential(scale=1/s)

    # advance time
    t = t + w

    # event occurs ()
    waiting_time.append(w)
    dinner_time.append(t)


# Convert to numpy arrays
waiting_time = np.array(waiting_time)
dinner_time = np.array(dinner_time)

# Event numbers for plotting
prey = np.arange(1, len(dinner_time) + 1)

print(waiting_time)
print(dinner_time)
print(len(dinner_time) - 1) #number prey found

plt.figure()
plt.step(dinner_time, prey, where='post')
plt.xlabel("Time")
plt.ylabel("Cumulative prey found")
plt.title("Individual searching: Poisson process")


# Poisson distribution
# simulate from above
