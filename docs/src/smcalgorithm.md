# SMC algorithm and particle approximations

## [SMC algorithm](@id basicsmc)

The implemented algorithm is as follows. The random variables $\zeta_p^1, \ldots, \zeta_p^N$ are typically referred to as particles, and the approximations they define are referred to as particle approximations.

---
1. Sample $\zeta_1^i \overset{\mathrm{i.i.d.}}{\sim} M_1$ for $i \in \{1,\ldots,N\}$.

1. For $p=2,\ldots,n$:

   * sample a vector $A_{p-1}^1,\ldots,A_{p-1}^N$ of i.i.d. ${\rm Categorical}(G_{p-1}(\zeta_{p-1}^1),\ldots,G_{p-1}(\zeta_{p-1}^N))$ random variables in increasing order.

   * sample $\zeta_p^i \overset{\mathrm{ind}} {\sim} M_p(\zeta_{p-1}^{A_{p-1}^i}, \cdot)$ for $i \in \{1,\ldots,N\}$.
---

This algorithm is referred to as SMC with multinomial resampling because the vector of numbers of *offspring* of the time $p-1$ particles follows a multinomial distribution.

## [Main particle approximations](@id mainapprox)

For some $p \in \{1,\ldots,n\}$, the particle approximation of $\eta_p(f)$ is obtained by constructing a particle approximation $\eta_p^N$ of $\eta_p$ and then calculating $\eta_p^N(f)$. Specifically, we define
```math
\eta_p^N := \frac{1}{N} \sum_{i=1}^N \delta_{\zeta_p^i}, \qquad p \in \{1,\ldots,n\}
```
and the particle approximation of $\eta_p(f)$ is then simply
```math
\eta_p^N(f) = \int f(x) \eta_p^N({\rm d}x) = \frac{1}{N} \sum_{i=1}^N f(\zeta_p^i).
```

Similarly, we define
```math
\hat{\eta}_p^N:= \frac{\eta_p^N \cdot G_p}{\eta_p^N(G_p)}, \qquad p \in \{1, \ldots, n\}
```
so that the particle approximation of $\hat{\eta}_p(f)$ is
```math
\hat{\eta}^N_p(f) = \int f(x) \hat{\eta}_p^N({\rm d}x) \frac{\sum_{i=1}^N G_p(\zeta_p^i) f(\zeta_p^i)}{\sum_{i=1}^N G_p(\zeta_p^i)}.
```

Finally, the particle approximations of $\hat{Z}_1,\ldots,\hat{Z}_n$ are defined by $\hat{Z}_1^N:=\eta_1^N(G_1)$ and for $p \in \{2,\ldots,n\}$, $\hat{Z}_p^N:=\hat{Z}_{p-1}^N\eta_p^N(G_p)$.

## Genealogical structure

The particle system defined by the SMC algorithm has a genealogical structure, induced by the ancestor indices. We define Eve indices such that $E_p^i$ is the index of the time $1$ ancestor of $\zeta_p^i$. That is, the Eve indices satisfy $E_1^i := i$ for $i \in \{1, \ldots, N\}$ and
```math
E_p^i = E_{p-1}^{A_{p-1}^i}, \qquad p \in \{2, \ldots, n\}, \quad i \in \{1, \ldots, N\}.
```
These can be used to quantify and visualize the path degeneracy phenomenon associated with the SMC algorithm. Since the ancestor indices are in sorted order, so too are the Eve indices at each time.

## Recursive definition

The SMC algorithm is intimately connected to the [recursive definition of the SMC measures](@ref recursivedefmeasures) themselves.

Let us define the sampling operator $S^N$ to be the random map that takes a measure $\mu$ and outputs a random measure
```math
\mu^N = S^N \mu = \frac{1}{N} \sum_{i=1}^N \delta_{X^i}, \qquad \text{where } X^i \overset{\mathrm{i.i.d.}}{\sim} \mu.
```
Then we have
```math
\eta^N_1 = S^N M_1, \qquad \hat{\eta}_1^N = \frac {\eta^N_1 \cdot G_1} {\eta^N_1(G_1)},
```
and for $p = 2, \ldots, n$,
```math
\eta^N_p = S^N (\hat{\eta}_{p-1}^N M_p), \qquad \hat{\eta}_p^N = \frac {\eta^N_p \cdot G_p} {\eta^N_p(G_p)}.
```

This mirrors almost exactly the recursive definitions for the SMC measures, where one has instead
```math
\eta_1 = M_1, \qquad \hat{\eta}_1 = \frac {\eta_1 \cdot G_1} {\eta_1(G_1)},
``` and for $p = 2, \ldots, n$,
```math
\eta_p = \hat{\eta}_{p-1} M_p, \qquad \hat{\eta}_p = \frac {\eta_p \cdot G_p} {\eta_p(G_p)}.
```

We notice also that $\hat{Z}^N_p = \hat{Z}^N_{p-1} \eta^N_p(G_p)$ while $\hat{Z}_p = \hat{Z}_{p-1} \eta_p(G_p)$ for $p \in \{2, \ldots, n\}$.
