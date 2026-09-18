# Continuous-time birth-death IBM

import numpy as np
import matplotlib.pyplot as plt
from time import time

# Function to grow the size of arrays
# Factor=2 (doubling size) by default
# Suggested fills for various dtypes
#     dtype    fill
#    ----------------------
#     float    0 or np.nan
#     int      -1
#     bool     False

def grow(x, fill, factor=2):

    # get original array size
    n = len(x)

    # make new bigger array (double by default)
    x_new = np.full(factor * n, fill, dtype=x.dtype)

    # copy original array into new one
    x_new[:n] = x

    return x_new



# Simulation control
rng = np.random.default_rng(95756)
t = 0.0      #needs to be float, e.g 0.0
t_max = 100


# Initial biological state

# --individual scale (parameters to set)
N_0 = 10
alive = np.full(N_0, True)  # alive state
b = np.full(N_0, 0.25)      # birth rate
m = np.full(N_0, 0.2)       # death rate

# --population scale (automatic)
N_t = np.sum(alive) #scale transition (sum over individuals -> population)
N = np.array([N_t])     #track N
times = np.array([t])   #track time

# --calculate initial event intensity
total_birth = np.sum(b * alive)
total_death = np.sum(m * alive)
intensity = total_birth + total_death

# Grow arrays to make room for new data

alive = grow(alive, fill=False, factor=100)
b = grow(b, fill=0, factor=100)
m = grow(m, fill=0, factor=100)
isize = len(alive)

N = grow(N, fill=-1, factor=1000)
times = grow(times, fill=np.nan, factor=1000)
psize = len(N)

# Simulate events

k = 0        #count number of events
irow = N_0   #array row we insert the next individual info into
prow = 1     #array row we insert the next population info into
chkpt_t = 1

while t < t_max and N_t > 0:

    # grow arrays if needed
    if irow == isize:
        alive = grow(alive, 0)
        b = grow(b, 0)
        m = grow(m, 0)
        isize = isize * 2

    if prow == psize:
        N = grow(N, -1)
        times = grow(times, np.nan)
        psize = psize * 2

    start_time = time() #timer

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
        # insert new individual in arrays
        alive[irow] = 1
        b[irow] = b[i]  #inherit from parent
        m[irow] = m[i]
        # update intensities
        total_birth = total_birth + b[i]
        total_death = total_death + m[i]
        # update N
        N_t = N_t + 1
        # increment row index for next time
        irow += 1
    else:
        # death
        alive[i] = 0
        # update intensities
        total_birth = total_birth - b[i]
        total_death = total_death - m[i]
        # update N
        N_t = N_t - 1

    intensity = total_birth + total_death

    # insert population level states in arrays
    N[prow] = N_t
    times[prow] = t
    # increment row index for next time
    prow += 1

    end_time = time() #stop timer

    k += 1 #event counter

    # some monitoring at integer time checkpoints
    if t > chkpt_t:
        tt = end_time - start_time
        print(f"N:{N_t} L:{len(alive)} t:{t:.1f} Evt:{k} Comp_t: {tt}")
        chkpt_t += 1


# Trim arrays (remove the unused elements)
irow = irow - 1
alive = alive[:irow]
b = b[:irow]
m = m[:irow]

prow = prow - 1
N = N[:prow]
times = times[:prow]

# Now we have clean arrays as before and can see results

print(alive)

plt.figure()
plt.step(times, N, where='post')
plt.ylim(bottom=0)
plt.xlabel("Time")
plt.ylabel("N")
plt.title("Birth-death IBM")
