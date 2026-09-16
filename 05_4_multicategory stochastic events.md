# Multicategory stochastic events


We have worked so far with binary stochastic events. However, many
biological processes involve more than two discrete options. For
example, does a bird lay a clutch of 0, 1, 2, or 3 eggs but never more?
Which of 5 nearest neighboring habitat patches does an individual
disperse to? Which of n possible partners does an individual mate with?
In situations like these the probabilities are usually unequal. For
example, 2 eggs could have the highest probability, dispersal
probability could depend on distance to the patch, mating probability
could depend on territory proximity and individual characteristics. The
different events are discrete and mutually exclusive and one event must
occur, so the probabilities necessarily sum to 1.

## Simulating a single multicategory stochastic event

The situations above are covered by the **categorical distribution**,
which is a generalization of the Bernoulli distribution from binary to
multiple categories. The categorical distribution is in turn a special
case of the multinomial distribution, where the number of trials
equals 1. To generate random draws from a categorical distribution we
use the multinomial distribution. For example, consider which of 5
nearest habitat patches an individual disperses to (5 event categories),
each with a different probability:

``` python
patch_ID = [35, 12, 47,  9, 37]
p = [0.1, 0.3, 0.2, 0.15, 0.25]
print(f"Total probability: {sum(p)}")
```

    Total probability: 1.0

To draw one stochastic dispersal event:

``` python
import numpy as np
rng = np.random.default_rng(16095)
outcome = rng.multinomial(n=1, pvals=p)
print(f"Categorical outcome: {outcome}")
```

    Categorical outcome: [0 0 0 1 0]

In this stochastic realization, the fourth category is the outcome and
thus the individual disperses to patch 9. The multinomial distribution
returns the number of events occurring in each category. Since we set
the number of trials equal to 1, one category will have a count of 1
while the others are 0. We can extract the index of the category using
`np.argmax()`:

``` python
event_index = np.argmax(outcome)
print(event_index)
```

    3

The event_index is an offset index, so here index 3 is the fourth event
category. The event (which patch dispersed to) can then be extracted:

``` python
disperse_to_patch = patch_ID[event_index]
print(f"Individual disperses to patch {disperse_to_patch}")
```

    Individual disperses to patch 9

While the above generates a categorical stochastic event from the
theoretical categorical distribution, an equivalent and more concise way
is to draw a random sample of size 1 from a list of the categories. This
is the most common approach.

``` python
disperse_to_patch = rng.choice(patch_ID, size=1, p=p)
print(f"Individual disperses to patch {disperse_to_patch}")
```

    Individual disperses to patch [12]

## Simulating multiple multicategory stochastic events

Just as the binomial distribution is the sum of multiple binary
stochastic events, the multinomial distribution is the sum of multiple
multicategory stochastic events. The multinomial distribution is useful
for simulating certain scenarios. For example, say we have 107
individuals to disperse to neighboring patches, then a stochastic
realization for the number of individuals dispersing to each patch can
be done with:

``` python
outcome = rng.multinomial(n=107, pvals=p)
print(f"Multinomial outcome: {outcome}")
print(f"Total dispersed: {sum(outcome)}")
print(f"Number dispersing to patch ID {patch_ID} is respectively {outcome}")
```

    Multinomial outcome: [17 26 24 17 23]
    Total dispersed: 107
    Number dispersing to patch ID [35, 12, 47, 9, 37] is respectively [17 26 24 17 23]
