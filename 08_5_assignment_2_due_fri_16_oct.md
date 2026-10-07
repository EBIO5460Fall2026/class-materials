# Assignment 2

**Due:** Friday 16 Oct 5:00 PM

**Grading criteria:** Complete the rubrics. On time submission.

**Percent of grade:** 7%

**Format for submitting assignments:**

* Submit code and answers to questions as comments in the same script
* You may use Python or R
* All code and text in one file please
* Please include your script (.py, .R) plus all the produced plots (.png, .jpg, .pdf). I don't want to run your code to make plots! Alternatively, you may produce a report (e.g. from a .rmd, .qmd, or other notebook) in one of the following formats: .md, .html, .pdf.
* The filename should be `pm4e_assignment2.R` or `pm4e_assignment2.py`. This will help me find the assignments in your repository.
* Supplemental files should be named similarly (e.g. `pm4e_assignment2_plot1.png`.

**Push your files to your GitHub repository**

## Learning goals

* Adapt the basic Gillespie individual-based algorithm to a new ecological scenario
* Reason about probability from individual-level stochastic events
* Design and conduct simulation studies to investigate ecological scenarios
* Connect individual-scale events to emergent group-level outcomes
* Recognize and reuse common modeling and algorithmic patterns across biological problems

## Building a many-eyes, co-operative foraging IBM

In this assignment you'll adapt the continuous-time birth-death IBM to model anti-predator vigilance in animal groups. The goal is not to write a new Gillespie algorithm from scratch. Instead, modify the birth-death IBM from week 5 Tuesday so that the biological processes and event types correspond to vigilance behavior. Don't optimize the model for speed. The goal is to first make a biologically intuitive prototype, even if it's slower.

### Ecological scenario

Consider a group of animals where individuals can either:

* Feed
* Remain vigilant for predators

Feeding allows an individual to acquire food, while vigilance allows an individual to detect approaching predators. Predators frequently attempt to prey on individuals in the group. If a vigilant individual detects a predator, it gives an alarm signal and all the individuals in the group are able to avoid predation. This is the "many eyes" hypothesis: a group may detect predators more effectively because multiple individuals contribute vigilance. How an individual allocates time feeding versus vigilance determines its fitness. As a simplified modeling exercise for learning, to minimize coding and computation, here we'll crudely define fitness as the multiple of survival and the amount of food an individual accumulates. If an individual succumbs to a predator, its fitness is zero. But there is a tradeoff. If an individual spends too much time being vigilant, it spends less time foraging. Maximum fitness is achieved by finding the right balance between vigilance and foraging. In reality, fitness is more rigorously defined than this, and individuals could adopt a variety of strategies, including cheating or cooperating. The simple model you develop here is easily extended to these more realistic considerations.

### Individual states

In the birth-death IBM, an individual's state was:

```python
alive = True or False
```

This is no longer sufficient. Add a second state variable that records whether an individual is:

```python
vigilant = True or False
```

If `vigilant` is `False`, the individual is feeding (unless it is no longer alive).

Track these group-level states through time:

* Number feeding, F
* Number vigilant, V

### Event types

The original birth-death IBM had two event types:

* Birth
* Death

Replace these with four new event types:

* Feeding event
* Predator encounter
* Switch feeding -> vigilant
* Switch vigilant -> feeding

You'll need to modify what happens during each event and update state variables. See below for ideas.

### Event intensities

In the birth-death IBM:

```python
birth intensity = total_birth
death intensity = total_death
```

For the vigilance model, derive expressions for:

* Feeding intensity
* Predator-encounter intensity
* Switch feeding -> vigilant intensity
* Switch vigilant -> feeding intensity

Hint: If feeding occurs at rate `f` per feeding individual and there are `F` feeding individuals, the total feeding intensity scales with `F`. Likewise, if individuals switch from feeding to vigilance at rate `switch_fv`, the total switching intensity will depend on the number currently feeding.

### Gillespie algorithm

The Gillespie algorithm is the same. After calculating event intensities:

1. Calculate the total intensity
2. Draw a waiting time from an exponential distribution
3. Determine which event type occurred
4. Determine which individual experienced the event (if required)
5. Update the system state

Only the number and types of events changes, as below.

### Feeding events

When a feeding event occurs:

* Choose one feeding individual
* Increase that individual's accumulated food by one unit

Add an array that tracks food accumulated by each individual.

### Behavior switching events

When a behavior switching event occurs:

* Choose an individual from among those in the appropriate current state
* Switch that individual's state

### Predator encounters

Predators attempt to attack the group at rate `a`. Suppose there are `V` vigilant individuals and each detects an attack independently with probability `p_detect`. Derive the probability `p_alarm` that the predator's attempt is detected by at least one vigilant individual. 

Hint: The logic is similar to the force-of-infection derivation in Assignment 1.

If the predator is not detected:

* Choose an individual from among those foraging (if any are foraging)
* That individual dies

### Update group-level states

Following each event:

* Recalculate and record `V`
* Recalculate and record `F`
* Record times of events

The variables `V` and `F` now play a role similar to `N` in the birth-death IBM.

### Model checking

Before proceeding to simulation studies, verify that your model behaves sensibly. Useful checks include:

* Are V + F and the number alive always equal?
* Do the event probabilities always sum to 1?
* Does increasing `p_detect` increase the frequency of successful alarms?
* Do feeding individuals accumulate more food than vigilant individuals?
* Plot V and F through time. Do the dynamics make sense?

### Encapsulate as a function

Now that you have a working model turn it into a function if you haven't already. Check the function works correctly.

It should return:

* Alive status for each individual
* Food accumulated for each individual
* Event times
* V
* F

### Rubric

- [ ] Individual state variables modified appropriately
- [ ] Four event types implemented
- [ ] Event intensities derived and implemented
- [ ] Feeding events implemented
- [ ] Behavior switching events implemented
- [ ] Predator encounter events implemented
- [ ] Derivation of p_alarm included; use text comments
- [ ] Model checked and debugged
- [ ] Model packaged as a function
- [ ] Function returns alive, food, times, V and F


## Simulation study 1

Simulate the model for 300 minutes (5 hours), corresponding to a day's worth of foraging activity. Use the following parameter values:

```python
N = 20            # group size
f = 1/2           # feeding rate, about 1 item every 2 minutes
a = 1/10          # predator arrival rate, about every 10 minutes
p_detect = 0.5    # probability a vigilant individual detects a predator
switch_fv = 1/5   # feeding -> vigilant, about every 5 minutes
switch_vf = 1     # vigilant -> feeding, about 1 minute in the vigilant state
```

For a single run of the model:

1. Plot the dynamics of number alive, V and F over time
2. Calculate the mean among times of the number of vigilant individuals. Compare to the expected number `switch_fv / (switch_fv + switch_vf)`. 
3. Calculate the mean fitness among individuals at the end of the realization. This is the per capita fitness and is equal to `mean(alive * food items accumulated)`.
4. Try several realizations and describe what happens

### Rubric

- [ ] Plot of number alive, V and F through time included (all on one plot)
- [ ] Mean number vigilant calculated
- [ ] Expected proportion vigilant calculated
- [ ] Simulated and expected vigilant compared
- [ ] Mean per-capita fitness calculated
- [ ] Behavior of the model described; use text comments


## Simulation study 2: fitness distribution

Simulate the model many times (e.g. 1000 realizations) using the same parameter values as in Simulation Study 1.

For each realization calculate:

```python
fitness = mean(alive * food)
```

where fitness represents per-capita food accumulated after accounting for mortality.

Store the fitness from every realization.

Hint: Python: appending to python base arrays is efficient to keep results. R: preallocate storage of the known size.

### Distribution of fitness

Plot the distribution of fitness among realizations.

Report:

* Mean fitness
* Standard deviation of fitness
* Monte Carlo standard error of the mean fitness

### Rubric

- [ ] Distribution of fitness plotted
- [ ] Mean fitness reported
- [ ] Standard deviation of fitness reported
- [ ] Monte Carlo standard error reported


### Simulation study 3: effect of vigilant behavior

The behavior-switching rates determine the long-run proportion of individuals that are vigilant and the amount of time an individual is vigilant vs feeding. Investigate how fitness changes as vigilance behavior changes.

Keep:

```python
switch_vf = 1
```

fixed and vary:

```python
switch_fv
```

across a grid of values from 1/5 to 1. Include at least 50 values in the grid.

For each value:

1. Run many realizations
2. Estimate mean fitness
3. Estimate standard error of fitness

Hints:

* You could use nested for loops, where the inner loop is the same as in simulation 1
* Time a small number of simulations to plan the tradeoff between number of simulations and time needed
* Use enough simulations to reduce the Monte Carlo error sufficiently to be confident in the conclusion

Plot:

* Mean fitness (and its 95% confidence interval) versus `switch_fv`

Code:

* Extract the value of `switch_fv` that maximized fitness in your simulation study

### Interpreting the many-eyes effect

Using comments in your script:

* Interpret the optimum value for `switch_fv`. About how often does an individual need to check for predators?
* Why does an intermediate level of vigilance maximize fitness?
* Is the vigilance level that maximizes individual fitness necessarily the same as the vigilance level that maximizes group survival?

### Rubric

- [ ] Simulation compares across multiple vigilance levels
- [ ] Plot of mean fitness versus `switch_fv` included
- [ ] 95% confidence interval for mean fitness included on the plot
- [ ] Biological interpretation of the vigilance-foraging tradeoff discussed