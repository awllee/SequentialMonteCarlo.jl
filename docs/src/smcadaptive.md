# [SMC with adaptive resampling](@id adaptiveresampling)

There is a tuning parameter associated with the ```SequentialMonteCarlo.smc``` algorithm that results in an adaptive version of the [SMC Algorithm](@ref basicsmc). This is the ```essThreshold``` parameter, whose use was proposed by Kong et al. (1994) and Liu & Chen (1995). We represent this parameter as $\tau$ below.

## (Relative) effective sample size

The adaptive SMC algorithm essentially involves choosing $(A_{p-1}^{1}, \ldots, A_{p-1}^N) = \left(1, \ldots, N \right)$ when the *weights* associated with particles $\zeta_{p-1}^{1}, \ldots, \zeta_{p-1}^N$ have a *relative effective sample size* exceeding $\tau$. This is sometimes referred to as "not resampling". The relative effective sample size of a collection of weights is defined as \\[ {\rm rESS}(w_{1}, \ldots, w_{N}) := \frac{\left(\frac{1}{N}\sum_{i=1}^N w_{i} \right)^{2}}{\frac{1}{N} \sum_{i=1}^N w_{i}^{2}}. \\] This function is invariant to rescaling of all of the weights by a common constant.

In the [SMC Algorithm](@ref basicsmc) we can view the weights associated with the particles $\zeta_{p-1}^{1}, \ldots, \zeta_{p-1}^N$ as being $W_{p-1}^{1}, \ldots, W_{p-1}^N$ where $W_{p-1}^{i} \propto G_{p-1}(\zeta_{p-1}^{i})$, and their relative effective sample size is ${\rm rESS}(W_{p-1}^{1}, \ldots, W_{p-1}^N)$.

In the SMC with adaptive resampling algorithm, the weight of a particle may be proportional to a product of many potential function values, depending on when the most recent normal resampling step was. In particular, whenever the standard resampling step is not taken, particles first inherit the weights of their ancestors and then multiply them by their own potential function values. The algorithm is as follows:

## Adaptive resampling algorithm

1. Sample $\zeta_{1}^{i} \overset{\mathrm{i.i.d.}}{\sim} M_{1}$ and compute $W_{1}^{i} \propto G_{1}(\zeta_{1}^{i})$ for $i\in \{1, \ldots, N\}$.

1. For $p=2,\ldots,n$:
    * compute $\mathcal{E}_{p-1}^N := {\rm rESS}(W_{p-1}^{1},\ldots,W_{p-1}^N)$.

    * if $p = n$ or $\mathcal{E}_{p-1}^N \leq \tau$ set $R_{p-1} \leftarrow 1$; otherwise set $R_{p-1} \leftarrow 0$.

    * if $R_{p-1} = 1$, sample a vector $A_{p-1}^{1}, \ldots, A_{p-1}^N$ of i.i.d. ${\rm Categorical}(W_{p-1}^{1}, \ldots, W_{p-1}^N)$ random variables in increasing order; otherwise set $(A_{p-1}^{1}, \ldots, A_{p-1}^N) = \left(1, \ldots, N \right)$.

    * sample $\zeta_{p}^{i} \overset{\mathrm{ind}}{\sim} M_{p}(\zeta_{p-1}^{A_{p-1}^{i}}, \cdot)$ and compute $W_{p}^{i} \propto \left(W_{p-1}^{i} \right)^{\mathbb{I}(R_{p-1} = 0)} G_{p}(\zeta_{p}^{i})$ for $i\in \{1,\ldots,N\}$.
---

Note that the time $n-1$ particles are always resampled, so that $\eta_{n}^N$ is always an unweighted particle approximation of $\eta_{n}$. This is primarily an implementation detail.

## Particle approximations

We define \\[ \eta_{p}^N \propto \sum_{i=1}^N \left(W_{p-1}^{i} \right)^{ \mathbb{I}(R_{p-1} = 0)} \delta_{\zeta_{p}^{i}}, \qquad p \in \\{1,\ldots,n\\}, \\] and \\[ \hat{\eta}_{p}^N \propto \sum_{i=1}^N W_{p}^{i} \delta_{\zeta_{p}^{i}}, \qquad p \in \\{1,\ldots,n\\}. \\]

Appropriate approximations of $\hat{Z}_{1}^N, \ldots, \hat{Z}_{n}^N$ are also well-defined, but tedious to display.

These particle approximations all enjoy the same [theoretical properties](@ref maintheory) stated for the standard SMC Algorithm, although the asymptotic variance maps $\sigma^2_p$ are generally different. In fact, adaptive resampling improves in certain scenarios the quality of SMC approximations; its effects have been studied theoretically by Del Moral et al. (2010) and Whiteley et al. (2016).

## Variance estimation

All of the methods described in [Variance estimators](@ref) can be run on output from the SMC with adaptive resampling algorithm.

One should be aware that the length of the vector returned by ```SequentialMonteCarlo.vpns``` will be of length $m = 1 + \sum_{i=1}^{n-1} R_i \leq n$. This is a consequence of the fact that resampling only at certain times can be viewed as running an SMC algorithm with modified Markov kernels and potential functions defined on a more complicated state space --- the details are not presented here --- so that the number of terms in the asymptotic variance decomposition is $m$ and not (necessarily) $n$.
