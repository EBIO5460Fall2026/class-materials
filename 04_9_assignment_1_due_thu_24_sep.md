# Assignment 1

**Due:** Thursday 24 Sep 3:30 PM

**Grading criteria:** Answer all the questions completely. On time submission.

**Percent of grade:** 7%

**Format for submitting assignments:**

* Submit code and answers to questions as comments in the same script
* You may use Python or R
* All code and text in one file please (i.e. not a separate file per question)
* Please include your script (.py, .R) plus all the produced plots (.png, .jpg, .pdf). I don't want to run your code to make plots! Alternatively, you may produce a report (e.g. from a .rmd, .qmd, or other notebook) in one of the following formats: .md, .html, .pdf.
* The filename should be `pm4e_assignment1.R` or `pm4e_assignment1.py`. This will help me find the assignments in your repository.
* Supplemental files should be named similarly (e.g. `pm4e_assignment1_plot1.png`.

**Push your files to your GitHub repository**


## Learning goals

* Adapt existing simulation code to new biological processes
* Practice reasoning and modeling using basic probability
* Connect individual-level stochastic events to population-level distributions
* Interpret probability distributions as emergent outcomes of biological mechanisms
* Use histograms from simulations to represent distributions
* Compare theoretical distributions to simulations

## Q1: IBM for binary stochastic events 

Modify the `death_IBM` model from `03_code_tue` to simulate a parasitic infection process instead

* Population size: N = 200
* Infected individuals: I = 20
* The remainder of the population is susceptible to infection: S = N - I
* Assume only individuals infected at the start of the year can transmit infection during that year. Newly infected individuals cannot infect others until the following year.
* Probability, $p_{ij}$, that infected individual $j$ infects individual $i$ during a one-year interval: 0.02
* Use a modern RNG for a Bernoulli distribution directly, instead of the hand coded one

Reasoning about probability

* You could simulate transmission of infection between each pair of individuals (but don't do that). Instead, we will scale up theoretically from pairwise transmissions to the "force of infection" from the population overall.
* Explain why the probability that a susceptible individual becomes infected is

$$
1 - (1 - p_{ij})^I
$$

* where $I$ is the number of infected individuals at the start of the year
* What assumptions are we making?

Investigations

* Track the random variables:
   * X = number of susceptible individuals that become infected after one year
   * Y = number of infected overall, including those initially infected
* Repeat simulations many times (e.g. 100,000) to estimate the distributions
* Plot the distribution for X (use as many simulations as needed for good resolution)
* Compare this simulation outcome to the appropriate theoretical distribution for X
* Compare simulated and theoretical means, variances, and standard deviation (sqrt(Var)) for X
* Plot the distribution for Y
* Compare simulated and theoretical means, variances, and standard deviation (sqrt(Var)) for Y

Tips and hints

* Begin by developing the model line by line before packaging as a function
* The distribution for Y, the number infected overall, is just the distribution for X offset by the number initially infected. Adjust the parameters of the distribution to account for the offset

Optional challenge

* Show that an IBM for pairwise infections is the same as the IBM with the scaled-up force of infection

Rubric

- [ ] Using text comments, explanation for force of infection included
- [ ] Using text comments, assumptions discussed
- [ ] Modified code included
- [ ] Bernoulli RNG used
- [ ] Plot included for distribution of X (susceptible to infected), with theoretical distribution overlaid
- [ ] Using text comments, simulated and theoretical mean, var, s.d. reported and compared for X
- [ ] Plot included for distribution of Y (total infected), with theoretical distribution overlaid
- [ ] Using text comments, simulated and theoretical mean, var, s.d. reported and compared for Y


## Q2: Reasoning about model assumptions

In the derivation of the exponential distribution we assumed

Q(t+Δt) = Q(t)Q(Δt)

Explain in your own words:

* What biological assumption justifies this step?
* What does "memoryless" mean in this context?
* Give an ecological or evolutionary example where this assumption might be reasonable
* Give an ecological or evolutionary example where it is unlikely to hold

Maximum 300 words


## Q3: Continuous time stochastic events

Part A: Modify the `death_continuous_time` model from `04_3_code_tue` to simulate an insect individual laying eggs

* Egg-laying rate: b = 0.1 per hour

Reasoning about probability

* What assumptions are we making about the egg laying process?

Investigations

* Simulate the outcome for one day, i.e. `t_max` = 24 hours
* Repeat simulations many times (e.g. 100,000) to estimate the distribution
* Plot the distribution for the number of eggs laid
* Compare this simulation outcome to the appropriate theoretical distribution
* Compare simulated and theoretical mean, var, s.d. for number of eggs laid

Part B: Now extend this model to incorporate a second step

* Egg laying events occur at the same rate as before, 0.1 events per hour
* Not all eggs are viable. An egg is viable with probability $p_v$ = 0.65

Investigations

* Simulate the outcome for one day, i.e. `t_max` = 24 hours
* Plot the distribution for the number of viable eggs laid
* Compare this simulation outcome to the appropriate theoretical distribution
* Compare simulated and theoretical mean, var, s.d. for number of viable eggs laid

Tips and hints

* The theoretical distribution in Part B is again $X \sim \text{Poisson}(\mu \times p_v)$, where $\mu$ is the mean of the original Poisson distribution. After Bernoulli "thinning", a Poisson process is still Poisson. This occurs in many situations, e.g.:
   * offspring born -> offspring surviving
   * animals present -> animals detected
   * infections present -> infections observed
   * seeds dispersed -> seeds germinating
   * mutation present -> function changed
* Use SciPy to calculate the theoretical PDF

Rubric

- [ ] Using text comments, assumptions discussed
- [ ] Modified code included
- [ ] Part A: Plot included for distribution of number of eggs laid, with theoretical distribution overlaid
- [ ] Part A: Using text comments, simulated and theoretical mean, var, s.d. reported and compared for eggs laid
- [ ] Part B: Plot included for distribution of number of viable eggs laid, with theoretical distribution overlaid
- [ ] Part B: Using text comments, simulated and theoretical mean, var, s.d. reported and compared for viable eggs laid


## Q4: New distribution - gamma from Poisson process

An individual needs to consume 5 prey items in a feeding bout to satisfy its physiological needs. Assuming the individual encounters prey at a rate of s = 0.1 per minute, how long does it take an individual to catch 5 prey items and complete a feeding bout? Modify the code in the `num_prey_found` model from `04_7_code_thu` so that the simulation stops when 5 prey items have been captured instead of stopping at a fixed time.

Investigations

* Repeat simulations many times (e.g. 100,000) to estimate the distribution
* Plot the distribution for length of feeding bout in minutes
* Compare the simulation outcome to the appropriate gamma distribution
* Compare simulated and theoretical mean, var, and s.d. for feeding bout length

Reasoning about probability

* What assumptions are we making about feeding bouts?
* How might ecological reality differ from these assumptions?

Tips and hints

* Consult wikipedia for how a gamma distribution arises, explanation of its parameters, means, and variances, and assumptions
* Use SciPy to calculate the theoretical PDF

Rubric

- [ ] Code included
- [ ] Plot included for distribution of feeding bout length, with theoretical distribution overlaid
- [ ] Using text comments, simulated and theoretical mean, var, s.d. reported and compared for length of feeding bout
- [ ] Using text comments, assumptions discussed

**Don't forget to use Piazza if you get stuck, to offer tips, share experience etc.**
