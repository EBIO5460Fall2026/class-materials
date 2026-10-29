# Assignment 3

**Due:** Friday 13 Nov 5:00 PM

**Grading criteria:** Complete the rubrics. On time submission.

**Percent of grade:** 11.67% (1/3 of assignment grade)

**Format for submitting assignments:**

* Submit code and answers to questions as comments in the same script
* You may use Python or R, or a mix if that's easier (R may be easier to follow my JAGS examples)
* All code and text in one file please (plus JAGS files), or two if you mix Py and R
* Please include your script (.py, .R, .jags) plus all produced plots (.png, .jpg, .pdf). I don't want to run your code to make plots! Alternatively, you may produce a report (e.g. from a .rmd, .qmd, or other notebook) in one of the following formats: .md, .html, .pdf.
* The filename should be `pm4e_assignment3.R` or `pm4e_assignment3.py`.
* Supplemental files should be named similarly (e.g. `pm4e_assignment3_plot1.png`).
* JAGS model files should be named sensibly (e.g. `bats_obs.jags`, `bats_ss.jags`)

**Push your files to your GitHub repository**

## Learning goals

* Develop a stochastic model from a biological and observation scenario
* Distinguish process and observation uncertainty
* Consider alternative models
* Understand the role of latent states in dynamic models
* Train a Bayesian dynamic model using MCMC via JAGS (training algorithm)
* Make model forecasts of dynamics including uncertainty (inference algorithm)
* Evaluate model assumptions and forecast performance

## Disease dynamics in a bat colony

In this assignment you'll develop and fit models for a pathogen outbreak in a bat colony. You might imagine a disease such as white nose syndrome. The pathogen causes infections, some bats recover while others die. Through the roosting season, ecologists conduct weekly surveys of the colony. Each week they:

* Estimate the size of the colony
* Sample bats and test them for infection

The data have already been collected and are provided in a CSV file.

This assignment is intended to bring together many of the ideas from the course. You'll now combine:

* Process models
* Observation models
* Bayesian training and inference
* Forecasting

## The data

The file is `bat_pathogen.csv` in the `data` directory.

The dataset contains weekly surveys with the following variables:

* Week
* Estimated colony size
* Number of infected bats detected in a sample
* Number of bats examined

We'll call these data the "real" data but they are simulated directly from the data generating process described below (thus, I know the underlying true process that you are trying to discover).

### First look at the data

Plot:

1. Colony-size estimate versus week
2. Observed prevalence versus week

where

prevalence = infected_sample / sample_size

* Does the population appear to be increasing or decreasing?
* Does infection prevalence appear to be increasing or decreasing?
* What biological processes might explain these patterns?

### Rubric

- [ ] Colony size plotted
- [ ] Infection prevalence plotted
- [ ] Patterns in colony size described
- [ ] Patterns in prevalence described
- [ ] Biological interpretation included


## Develop the disease model

Turn the following mathematical description of the biological processes into code.

We assume the colony is closed, i.e. no immigration or emigration. The colony migrates and roosts together. Individuals can be in one of four states:

* S, Susceptible
* I, Infected
* R, Recovered
* D, Dead

An individual can transition from state S to I, or from state I to either R or D. This is a SIRD model.

### Force of infection

Use the force of infection model from Assignment 1. The parameter

$\alpha$

is the transmission probability for a single infected bat. The force of infection is then:

$FOI_t = 1 - (1 - \alpha)^{I_t}$

where $I_t$ is the number infected.

### Infection

New infections ("**C**ases") are

$C_t \sim \text{Binomial}(S_t, FOI_t)$

### Infection resolution

The parameter

$q$

is the weekly probability that an infection resolves (to either recovered or dead).

Then

$L_t \sim \text{Binomial}(I_t, q)$

where $L_t$ is the number of infections that resolve ("**L**eave" the infected state) during a week.


### Recovery and mortality

The parameter

$m$

is the probability that a resolved infection ends in death.

Then deaths ("**M**ortality") are

$M_t \sim \text{Binomial}(L_t, m)$

and recoveries ("**G**ains" to the recovered state) are

$G_t = L_t - M_t$

### State updates

$S_{t+1} = S_t - C_t$

$I_{t+1} = I_t + C_t - L_t$

$R_{t+1} = R_t + G_t$

$D_{t+1} = D_t + M_t$

$N_t = S_t + I_t + R_t$

### Interpret the parameters

Describe in comments:

* Biological meaning of $\alpha$
* Biological meaning of $q$
* Biological meaning of $m$

### Rubric

- [ ] Infection, resolution, and mortality processes coded
- [ ] State updates coded
- [ ] Interpretation of three parameters included


## Observation model

The survey data does not directly observe the true state variables. States are measured with estimation or sampling error. Turn the following mathematical description of the observation processes into code.


### Infection survey

A sample of bats is examined (ideally a random sample). There is no error in detecting the pathogen if a bat is infected. The number infected is

$y_t \sim \text{Binomial} \left( n, \frac{I_t}{N_t} \right)$

where $n$ is the sample size and $I_t / N_t$ is the prevalence of the pathogen.


### Colony-size estimate

The colony size is not observed perfectly. The ecologists try to count all of them but bats move around and people lose count. We'll assume errors are proportional to the number of bats with a lognormal observation process (simulated data are rounded to the nearest integer):

$\hat N_t \sim \text{Lognormal} \left( \ln(N_t), \sigma_N \right)$

where $\sigma_N$ is the observation error.

### Rubric

- [ ] Infection observation process coded
- [ ] Colony-size observation process coded


# Simulation study 1: data generating process

Simulate the model for the full data generating process including process uncertainty (demographic stochasticity) and observation uncertainty.

### Explore the model

Simulate the model and plot:

True states:

* Susceptible bats
* Infected bats
* Recovered bats
* Cumulative deaths
* Colony size

Generated data:

* Observed prevalence
* Observed colony size

Experiment with parameter values.

Questions:

1. How does increasing transmission, $\alpha$, affect epidemic dynamics?
2. How does increasing mortality, $m$, affect colony size?
3. How does increasing the resolve probability, $q$, affect epidemic duration?
4. Do the simulated patterns qualitatively resemble the provided "real" data?

### Rubric

- [ ] Plot code included (no need to include plots)
- [ ] Parameters explored and described for Q1-3
- [ ] Qualitative descriptive comparison to data included


# Simulation study 2: observation-error-only model

Now code an analogous model that ignores demographic stochasticity. Replace the stochastic process terms with their expectations, e.g.

$C_t \sim \text{Binomial}(S_t,FOI_t)$

becomes

$C_t = FOI_t S_t$

The process model is now deterministic. Keep the same stochastic observation model.


## Compare the models

Using the same parameter values:

1. Simulate the stochastic process model
2. Simulate the deterministic process model
3. Compare the trajectories

Questions:

* Which model is more variable?
* Do they have similar average behavior?
* How much variability is generated by demographic stochasticity?

### Rubric

- [ ] Observation-error-only model implemented
- [ ] Stochastic and deterministic process models compared in text comments


# Bayesian inference using JAGS

Train the **observation-error-only model** using JAGS. Stan is more difficult for the full state-space model because of the discrete latent states.

Hint: the JAGS model code will closely resemble the simulation code.

## Suggested priors

Use:

```jags
alpha ~ dunif(0, 0.01)
q ~ dbeta(1, 1)
m ~ dbeta(1, 1)
sigma_N ~ dexp(1)
tau_N <- 1 / pow(sigma_N, 2)
lnN1 ~ dnorm(log(<N1_obs>), 0.3)
N1 <- exp(lnN1)
I1 ~ dnorm(0, 10) T(0,)
```
where <N1_obs> is the literal number of the colony size recorded for the first weekly visit.

Assume

```jags
R[1] <- 0
D[1] <- 0
```

## Check your model code first

Before fitting the real data:

1. Simulate a dataset from your original observation-error only model
2. Fit the model to the simulated data
3. Compare estimated parameters to the known values
4. Try this a few times

If the estimated parameters are about right, this provides evidence that the JAGS model behaves as intended.

### Rubric

- [ ] JAGS model coded
- [ ] Results discussed in comment text


## Train the model with the real data

Use the standard workflow:

1. Produce MCMC samples
2. Examine trace plots, Rhats, effective sample sizes
3. Adjust burn-in and iterations as needed
4. Perform posterior predictive checks

Report:

* Posterior mean
* 95% credible interval (CPI)

for $\alpha$, $q$, $m$, $\sigma_N$, $N[1]$, $I[1]$.

Hint: print the fit object

### Rubric

- [ ] Model fitted to data
- [ ] Results for trace plots, diagnostics, ppchecks discussed
- [ ] Posterior summaries reported as text comments


## Forecast the remainder of the season

We are at week 20 right now. The colony will eventually leave the cave and migrate but until then the rest of the roosting season is ahead. A forecast is a prediction about the future. How well can we forecast the rest of the season using this trained model?


### Generate forecast samples

Using posterior samples from the fitted model generate predicted trajectories for weeks 1 - 40. For each posterior sample, simulate these states of the biological system:

* Number infected
* Colony size

Also simulate:

* Observed number infected
* Observed colony size

You now have posterior samples of forecasts for trajectories of system states and for new data.

Hints:

* Reuse your simulation code for the observation-error only model (packaged into a function)
* Plug values for each posterior sample into the model


## Visualization

Make two plots: one for number infected, the other for colony size.

### States

Start by plotting states with their intervals. 

Number infected vs time

* Posterior mean trajectory
* 95% credible interval

Colony size vs time

* Posterior mean trajectory
* 95% credible interval

### Prediction intervals

A prediction interval shows what we would expect for collected new data. Add prediction intervals to the plot. Compute these from the posterior samples for "observed ...".

Observed number infected vs time

* 95% credible interval only (prediction interval)

Observed colony size vs time

* 95% credible interval only (prediction interval)

### Example trajectories

A useful visualization is to overlay 10–20 posterior trajectories of the states (randomly sample posterior draws because nearby MCMC samples are correlated).

### Data

Overlay the real data (compute $I$ from colony size and prevalence).

### Post forecast plots to Piazza for comparison

After posting your two plots, send me an email and I will send the remaining 20 weeks of data.

### Questions

* Does the fitted model describe the observed data?
* Does the forecast appear biologically reasonable?

## Rubric

- [ ] Forecast code included
- [ ] Number infected forecast plot submitted
- [ ] Colony size forecast plot submitted
- [ ] Credible intervals shown
- [ ] Two questions discussed in comment text
- [ ] Forecast posted to Piazza


## Forecast reveal and model update

After receiving the full dataset of 40 weeks:

1. Compare forecast versus observations
2. Retrain the model using all 40 weeks
3. Produce updated parameter estimates and trajectories

Questions:

* How accurate was the forecast compared to data?
* Which features of the epidemic were predicted well?
* Which features were predicted poorly?
* What assumptions of the model may have contributed to forecasting errors?

### Rubric

- [ ] Plot of forecast compared to observations submitted
- [ ] Model retrained with full dataset
- [ ] Updated plots submitted
- [ ] Forecast accuracy and limitations discussed as text comment



# Optional challenges


## Challenge 1: demographic stochasticity in forecast

Keep the fitted parameter values from the observation-error-only model. Now include demographic stochasticity in the forecast. No new parameters are required.

Questions:

* How does forecast uncertainty change?
* Does uncertainty increase through time?
* Why?


## Challenge 2: full state-space model

Write the JAGS model that includes both:

* Demographic stochasticity
* Observation error

Try training the model

* Use the full dataset.
* Start by fixing the initial values of N and I to values obtained from the observation-error-only fit.
* Use the same priors for the remaining parameters.
* Use `jags.parallel()` to speed computation.

Investigate:

* Trace plots
* Autocorrelation
* Effective sample size
* Convergence diagnostics
* Adjust burn in, sampling and thinning

Encouragement: the model will converge with sufficient burn in and sampling but it may take several hours.

Questions:

* How does the computation compare to the observation-error-only model?
* Why might the full state-space model be substantially harder to fit?
