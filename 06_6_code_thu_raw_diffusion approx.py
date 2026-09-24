import numpy as np
import matplotlib.pyplot as plt

b = 0.25
m = 0.2
N_0 = 10
dt = 0.01
t_max = 20

rng = np.random.default_rng()
N = [N_0]
times = [0]

t = 0
N_t = N_0
while t < t_max:
    t = t + dt
    Z = rng.normal()
    N_t = N_t + (b-m) * N_t * dt + np.sqrt((b+m) * N_t * dt) * Z
    times.append(t)
    N.append(N_t)

times = np.array(times)
N = np.array(N)

plt.figure()
plt.plot(times, N)
plt.ylim(bottom=0)
plt.xlabel("Time")
plt.ylabel("N")
plt.title("Birth-death stochastic approximation")


# Package into a function

def birth_death_diff_approx(N_0, b, m, dt, t_max, rng):

    N = [N_0]
    times = [0]

    t = 0
    N_t = N_0
    while t < t_max:
        t = t + dt
        Z = rng.normal()
        N_t = N_t + (b-m) * N_t * dt + np.sqrt((b+m) * N_t * dt) * Z
        if N_t < 0:
            N_t = 0
        times.append(t)
        N.append(N_t)

    times = np.array(times)
    N = np.array(N)

    return N, times


# Simulation study

b = 0.25
m = 0.2
N_0 = 10
dt = 0.01
t_max = 100

rng = np.random.default_rng()
N = []
for i in range(25):
    N_series, t = birth_death_diff_approx(N_0, b, m, dt, t_max, rng)
    N.append(N_series)


plt.figure()
for i in range(25):
    plt.plot(t, N[i], lw=0.75, alpha=0.5)
plt.ylim(bottom=0)
# plt.ylim(top=10) # to see extinctions
plt.xlabel("Time")
plt.ylabel("N")
plt.title("Birth-death stochastic approximation")
