# [Theoretical properties of particle approximations](@id maintheory)

All convergence results below involve taking limits as $N\rightarrow\infty$. Convergence almost surely is denoted $\overset{a.s.}{\rightarrow}$, convergence in probability is denoted $\overset{P}{\rightarrow}$ and weak convergence or convergence in distribution/law is denoted $\overset{L}{\rightarrow}$.

The following Theorem provides some justification for the use of the [main particle approximations](@ref mainapprox); these now classical results can be deduced from various results of Del Moral (2004). We present the results only for the "hatted" or "updated" quantities to avoid repetition; analogous hold for their "unhatted" counterparts with a different sequence of maps $\sigma_1^2, \ldots, \sigma_n^2$.

---
**Theorem** [Del Moral, 2004]. Let the potential functions $G_1, \ldots, G_n$ be bounded and strictly positive. There exist maps $\hat{\sigma}_1^2, \ldots, \hat{\sigma}_n^2$ such that the following hold for an arbitrary, bounded $f$:

* Lack-of-bias: $\mathbb{E} \left[ \hat{\gamma}_p^N(f) \right] = \hat{\gamma}_p(f)$ for all $N \geq 1$.

* Consistency: $\hat{\gamma}_p^N(f)\overset{a.s.}{\rightarrow}\hat{\gamma}_p(f)$ and $\hat{\eta}_p^N(f)\overset{a.s.}{\rightarrow}\hat{\eta}_p(f)$.

* Asymptotic variance and mean-squared error (MSE):

```math
N {\rm var} \left \{ \hat{\gamma}_p^N(f) / \hat{\gamma}_p(1) \right \} \rightarrow \hat{\sigma}_p^2(f),
```
  and
```math
N \mathbb{E} \left[ \left\{ \hat{\eta}_p^N(f) - \hat{\eta}_p(f) \right\} ^2 \right] \rightarrow \hat{\sigma}_p^2(f-\hat{\eta}_p(f)).
```

* Central Limit Theorems:
```math
\sqrt{N} \left( \hat{\gamma}_p^N(f) / \hat{\gamma}_p(1) - \eta_p^N(f) \right) \overset{L}{\rightarrow} N(0,\hat{\sigma}_p^2(f)),
```
  and
```math
\sqrt{N} \left( \hat{\eta}_p^N(f) - \eta_p^N(f) \right) \overset{L}{\rightarrow} N(0,\hat{\sigma}_p^2(f - \hat{\eta}_p(f))).
```
---

### Note on the sorted ancestor indices

The theoretical results above are typically proven for an algorithm that differs very slightly from the [SMC algorithm](@ref basicsmc) implemented here. In particular, one would usually analyze the algorithm by considering $A_{p-1}^1,\ldots,A_{p-1}^N$ to be i.i.d. ${\rm Categorical} (G_{p-1}(\zeta_{p-1}^1), \ldots, G_{p-1}(\zeta_{p-1}^N))$ random variables rather than being in sorted order.

The laws of the approximations $\hat{Z}_p^N$, $\eta_p^N(f)$ and $\hat{\eta}_p^N(f)$, however, are invariant to permutations of the indices of the ancestors $A_{p-1}^1,\ldots,A_{p-1}^N$ in the algorithm. Therefore, the sorting of the ancestor indices may be regarded as an implementation issue that does not affect the particle approximations themselves.
