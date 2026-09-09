# Deriving the exponential distribution from the algorithm

Assumptions

* small time interval $\Delta t$
* probability of death $m \Delta t$ during the interval
* constant probability through time

The parameter $m$ can be thought of as the mortality rate. As we'll see, its inverse $1/m$ is the expected age of an individual.

The probability of death of an individual during the interval $\Delta t$, given it is alive so far, is

$$
\text{Pr}(\text{death in } \Delta t \ | \ \text{alive so far}) = m \Delta t + o(\Delta t) \ \ \ \ \ \ \ \ \ \ \ \ \text{(Eq. 1)}
$$

where

$$
o(\Delta t) = a \Delta t^2 + b \Delta t^3 + \ ...
$$

The role of $o(\Delta t)$, pronounced "oh delta tee", is pretty neat. The idea is that we can use a **polynomial function** to represent the (potentially nonlinear) relationship between probability of death and interval duration to infinite accuracy. In Eq. 1 we are breaking out the linear term, $m \Delta t$, with $o(\Delta t)$ containing all the higher order polynomial terms. We can think of $m \Delta t$ as the **linear approximation**, or "first-order" approximation, of the probability and $o(\Delta t)$ as the error. In our simulation, we used only the linear approximation term.

As with many problems in probability, the **probability of an event not happening**, together with thinking somewhat backwards, is the key to progress. We often use the notation $q$ for the probability of an event not happening. The probability of no death in the interval, $Q(\Delta t)$, is one minus the probability of death

$$
\begin{align}
    Q(\Delta t) &= \text{Pr}(\text{no death in } \Delta t) \\
        &= 1 - m \Delta t + o(\Delta t). \ \ \ \ \ \ \ \ \ \ \ \ \text{(Eq. 2)}
\end{align}
$$

Now, to look ahead, at this point we're trying to make Eq. 4, so we need $Q(t)$ and $Q(t + \Delta t)$. To get there we ask: what is the probability of **not dying** from 0 to $t + \Delta t$?

First the individual must not die from 0 to $t$, which has probability

$$
Q(t) = \text{Pr}(\text{no death from } 0 \text{ to } t). \ \ \ \ \ \ \ \ \ \ \ \ \text{(Eq. 3)}
$$

Then it must also not die in the interval $\Delta t$ right after that, which has probability $Q(\Delta t)$ (Eq. 2). Assuming these two events are **independent**, which means they don't depend on time or each other (often said to be memoryless), the probability of not dying in either period is the **product** of their probabilities, that is

$$
\begin{align}
    Q(t + \Delta t) &= Q(t)Q(\Delta t) \\
        &= Q(t)(1 - m \Delta t + o(\Delta t)) \\
        &= Q(t) - m Q(t) \Delta t + o(\Delta t).
\end{align}
$$

Now comes a neat trick from calculus. Subtract $Q(t)$ from both sides and divide by $\Delta t$:

$$
\displaystyle \frac{Q(t + \Delta t) - Q(t)}{\Delta t} = \frac{-m Q(t) \Delta t}{\Delta t} + \frac{o(\Delta t)}{\Delta t} \ \ \ \ \ \ \ \ \ \ \ \ \text{(Eq. 4)}
$$

We can turn this equation into **continuous time** by shrinking $\Delta t$ to zero, which mathematically is taking the **limit** $\Delta t \rightarrow 0$. The left hand side of Eq. 4 together with the limit is the **definition of the derivative** of $Q(t)$ with respect to $t$. When we shrink $\Delta t$ to zero, the higher order terms on the right hand side vanish, and we get (by definition) the **differential equation**

$$
\displaystyle \frac{dQ(t)}{dt} = -m Q(t).
$$

If we then solve this differential equation (by integrating over $t$) to find $Q(t)$, we get

$$
Q(t) = Q(0)e^{-mt}.
$$

Now $Q(0)$ = 1 because we assume at time 0 the individual was alive. Thus we get

$$
Q(t) = e^{-mt}. \ \ \ \ \ \ \ \ \ \ \ \ \text{(Eq. 5)}
$$

This is an exponential equation, so we're starting to see a hint of why the distribution is exponential.

But we're not done yet, we must keep thinking backwards! The probability of no death up until time $t$ is the same as the probability that the time of death, $T_\text{death}$, is greater than $t$, that is

$$
Q(t) = \text{Pr}(T_\text{death} > t).
$$

Therefore, by the law of total probability, and substituting Eq. 5


$$
\begin{align}
    \text{Pr}(T_\text{death} \leq t) &= 1 - Q(t) \\
        &= 1 - e^{-mt} \ \ \ \ \ \ \ \ \ \ \ \ \text{(Eq. 6)}
\end{align}
$$

This equation is the cumulative distribution function (CDF) for the time to death. The CDF represents the total probability from $0$ to $t$, which is the **area under the curve** of the probability distribution function (PDF), or the **integral** of the PDF from $0$ to $t$. Conversely, the **derivative** of the CDF is the PDF, so taking the derivative of Eq. 6 gives us the PDF:

$$
f(t) = me^{-mt}. \ \ \ \ \ \ \ \ \ \ \ \ \text{(Eq. 7)}
$$

We have arrived! This is the PDF of an **exponential distribution**. The mean of this distribution, i.e the mean time to death (or age of an individual), is $1/m$. 

So, our discrete-time algorithm for time to death gives an exponential distribution (Eq. 7) as the duration of the time increment $\Delta t$ approaches zero. This matches our experience with simulating the algorithm.
