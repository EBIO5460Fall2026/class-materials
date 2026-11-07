# Discussion papers for week 13

You need to read the paper plus submit a small written reflection (except for the discussion leader). The reflection should be **pushed to github before class**.  You won't get credit for the discussion unless you push to github before class. The reflection should include your thoughts on each of the questions under "questions to consider" below. It can be very brief and in bullet-point form but feel free to expand if you prefer.

## Tuesday

* Leader: MM
* Kéry et al (2013). Analysing and mapping species range dynamics using occupancy models. **Journal of Biogeography** 40:1463-1474. https://doi.org/10.1111/jbi.12087.

## Thursday

* Leader: VP
* Buckley & Roughgarden (2006). Climate, competition, and the coexistence of island lizards. **Functional Ecology** 20:315-322. https://doi.org/10.1111/j.1365-2435.2006.01095.x.

## Questions to consider for the papers

Depending on the paper, not all questions are necessarily relevant or might need to be adapted, e.g. if there are multiple models or case studies, you could focus on one.

* What is the biological or ecological question being addressed?
  * How is a model helpful to answer the question?
* Classify the process model
  * broad: ordinary or partial differential equations, IBM, difference equations, etc
  * "marketing" name: state-space, IPM, DOM ...
  * stochastic/deterministic
  * discrete/continuous space/time
* Describe the assumptions of the model.
  * It can be helpful here to draw a conceptual model diagram (no need to submit but good idea for leaders to include in presentation)
  * What are the state variables? (quantities the model tracks)
  * What biological processes are included?
  * How are biological processes represented? (functions, algorithms, probability distributions)
  * If modeled above the individual level (e.g. as an ODE), how do functions relate to **events** that happen to **individuals** (or lower level units)?
  * How are observation processes represented?
  * If stochastic, what are the sources of stochasticity? (demographic, environmental, observation error, other)
* Identify and describe the parameters of the model
  * What biological quantity or process does each represent?
  * Which parameters are estimated, fixed, or obtained from previous studies?
* How are data connected to the model?
  * What data were collected?
  * What variables are directly observed versus latent or unobserved?
  * Is there an observation model?
  * How are the parameters estimated? (likelihood, Bayesian, SBI, other)
* What techniques are used to evaluate the model?
  * predictive performance, posterior predictive checks, model comparison etc
  * sensitivity analysis, stablity analysis etc
* What biological insights come from the model?
  * Which conclusions depend critically on the model?
  * Did the model reveal anything that would not have been obvious from the data alone?
* What is the scope of inference?
  * What population, system, species, place, or time period are the authors trying to generalize to?
  * What assumptions are required for that generalization?
* Most papers have some problems! Do you have any critiques of the paper, either methodological or about the conclusions drawn? e.g.
  * Are potentially important mechanisms omitted?
  * How realistic are assumptions?
  * Are there weaknesses/strengths in the chosen level of complexity?
  * Are there alternative models or explanations that should have been considered?
* 2-3 other bullet points, which might be insights or questions you want to raise in the discussion

## Discussion leaders

Everybody will consider the questions above for each paper. The discussion leader will in addition prepare a presentation (e.g. powerpoint, google slides, beamer) that answers each of the **first 8** questions (the what, how, and why). The leader should go into a little bit of detail for each question (e.g. one slide per question; no more than two slides) highlighting especially any method details we haven't already covered in the first part of the course.  The overview presentation will be about 15 minutes. The discussion leader should also think of 2-3 key questions to ask of the group as discussion points. These questions should go beyond summarizing the paper and instead encourage evaluation of the modeling choices, assumptions, inference methods, and biological conclusions. The discussion leader doesn't need to submit a written reflection but should submit their presentation.
