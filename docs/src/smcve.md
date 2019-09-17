# [Variance estimators](@id vepage)

## Asymptotic variances

We recall that the particle approximations have [asymptotic variances](@ref maintheory) associated with particular test functions $f$. Letting $n$ be arbitrary to avoid additional indexing notation, the values $\sigma_n^2(f)$ and $\hat{\sigma}_n^2(f)$ are typically not possible to calculate exactly. They can be decomposed into sums as follows:
```math
\sigma_n^2(f) = \sum_{p=1}^n v_{p,n}(f), \qquad \hat{\sigma}_n^2(f) = \sum_{p=1}^n \hat{v}_{p,n}(f).
```
The quantities $v_{p,n}(f)$ and $\hat{v}_{p,n}(f)$ also cannot be computed exactly in general, but this decomposition can in some circumstances shed some light on the nature of the approximation errors.

## The estimators

Most Monte Carlo approximations are accompanied by Monte Carlo approximations of their asymptotic variance or MSE. Lee and Whiteley (2015) defined $V_p^N(f)$, computable using the simulated particle system, allowing one to approximate asymptotic variances of interest.

| approximation | approximation of |
| ------------- | ---------------- |
| $V_p^N(f)$ | ``{\rm var} \left \{ \gamma_p^N(f)/\gamma_p(1) \right \}`` |
| $V_p^N(f-\eta_p^N(f))$ | ``\mathbb{E} \left[ \left\{ \eta_p^N(f) - \eta_p(f) \right\}^2 \right]`` |
| $\hat{V}_p^N(f)$ | ``{\rm var} \left\{ \hat{\gamma}_p^N(f) / \hat{\gamma}_p(1) \right\} `` |
| $\hat{V}_p^N(f-\hat{\eta}_p^N(f))$ | ``\mathbb{E} \left[ \left\{ \hat{\eta}_p^N(f) - \hat{\eta}_p(f) \right\} ^2\right]`` |

These approximations have the following theoretical justification, particularly for their approximation of the asymptotic variance and MSE maps defined earlier. Again, we present the results only for the "updated" quantities to avoid repetition; analogous results hold for their "unhatted" counterparts.

## Theoretical justification

The variance approximations have consistency properties, and some of them have lack-of-bias properties.

---
**Theorem** [Lee & Whiteley, 2015]. Let the potential functions $G_1, \ldots, G_n$ be bounded and strictly positive. The following hold for an arbitrary, bounded $f$:

1. Lack-of-bias: $\mathbb{E}\left[(\hat{Z}_p^N)^2\hat{V}_p^N(f)\right]={\rm var}\left\{ \hat{\gamma}_p(f)\right\}$ for all $N \geq 1$.

1. Consistency: $N\hat{V}_p^N(f) \overset{P}{\rightarrow} \hat{\sigma}_p^2(f)$ and $N \hat{V}_p^N(f-\hat{\eta}_p^N(f)) \overset{P}{\rightarrow} \hat{\sigma}_p^2(f - \hat{\eta}_p(f))$.
---

**Remark**. A consistent estimator of the asymptotic MSE of $\eta_p(f)$ or $\hat{\eta}_p(f)$ was proposed in Chan & Lai (2013). This is very similar to that of Lee & Whiteley (2015) but does not satisfy the lack-of-bias property.

Lee & Whiteley (2015) also define approximations $v_{p,n}^N(f)$ and $\hat{v}_{p,n}^N(f)$ of $v_{p,n}(f)$ and $\hat{v}_{p,n}(f)$, respectively. Their sums $v_n^N(f)=\sum_{p=1}^nv_{p,n}^N(f)$ and $\hat{v}_n^N(f)=\sum_{p=1}^n\hat{v}_{p,n}^N(f)$ can also be used as alternative approximations of $\sigma_n^2(f)$ and $\hat{\sigma}_n^2(f)$

---
**Theorem** [Lee & Whiteley, 2015]. Let the potential functions $G_1, \ldots, G_n$ be bounded and strictly positive. The following hold for an arbitrary, bounded $f$:

1. Lack-of-bias: $\mathbb{E}\left[\left(\hat{Z}_n^N\right)^2\hat{v}_{p,n}^N(f)\right]=\hat{Z}_n^2\hat{v}_{p,n}(f)$ for all $N\geq1$.

2. Consistency: $\hat{v}_{p,n}^N(f)\overset{P}{\rightarrow}\hat{v}_{p,n}(f)$ and $\hat{v}_{p,n}^N(f)(f-\hat{\eta}_n^N(f))\overset{P}{\rightarrow}\hat{v}_{p,n}(f-\hat{\eta}_n(f))$.

3. Lack-of-bias: $\mathbb{E}\left[\left(\hat{Z}_n^N\right)^2\hat{v}_n^N(f)\right]=\hat{Z}_n^2\hat{\sigma}_n^2(f)$ for all $N\geq1$.

4. Consistency: $\hat{v}_n^N(f)\overset{P}{\rightarrow}\hat{\sigma}_n^2(f)$ and $\hat{v}_n^N(f-\hat{\eta}_n^N(f))\overset{P}{\rightarrow}\hat{\sigma}_n^2(f-\hat{\eta}_n(f))$.
---
