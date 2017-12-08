# Conditional SMC

One relatively recent subtle modification of the SMC algorithm was proposed by Andrieu et al. (2010), and known as conditional SMC (cSMC). The algorithm is as follows:

---
Input: reference path $(x_{1}^{{\rm ref}},\ldots,x_{n}^{{\rm ref}})$

1. Set $\zeta_{1}^{1}=x_{1}^{{\rm ref}}$ and sample $\zeta_{1}^{i}\overset{\mathrm{i.i.d.}}{\sim}M_{1}$ for $i \in \{2, \ldots, N \}$.

1. For $p=2,\ldots,n$:

    * set $A_{p-1}^{1}=1$ and $\zeta_{p}^{1}=x_{p}^{{\rm ref}}$.

    * sample a vector $A_{p-1}^{2},\ldots,A_{p-1}^{N}$ of i.i.d. ${\rm Categorical}(G_{p-1}(\zeta_{p-1}^{1}),\ldots,G_{p-1}(\zeta_{p-1}^{N}))$ random variables in increasing order.

    * sample $\zeta_{p}^{i} \overset{\mathrm{ind}}{\sim} M_{p}(\zeta_{p-1}^{A_{p-1}^{i}}, \cdot)$ for $i \in \{2, \ldots, N \}$.

1. Sample $K_{n} \sim {\rm Categorical}(G_{n}(\zeta_{n}^{1}),\ldots,G_{n}(\zeta_{n}^{N}))$ and for $p=n-1,\ldots,1$ set $K_{p}=A_{p}^{K_{p+1}}$.

1. Output the path $(\zeta_{1}^{K_{1}},\ldots,\zeta_{n}^{K_{n}})$.
---

cSMC is of special interest as it defines a Markov kernel that is ergodic, for $N \geq 2$, with invariant probability measure given by \\[ \hat{\pi}(A)=\hat{Z}_{n}^{-1}\int_{A}G_{n}(x_{n})M_{1}({\rm d}x_{1})\prod_{p=2}^{n}G_{p-1}(x_{p-1})M_{p}(x_{p-1},{\rm d}x_{p}), \qquad A \in \mathcal{X}^{\otimes{n}}. \\]

That is, the Markov kernel $P_{N}(x^{{\rm ref}}, \cdot)$ defined by calling the cSMC algorithm has invariant probability measure $\hat{\pi}$. Theoretical results, in particular concerning the effect of $N$ on the rate of convergence of $P_{N}$, have been developed by Chopin & Singh (2015), Andrieu et al. (2018) and Lindsten et al. (2015).

Conditional SMC may also be used with adaptive resampling as detailed in [Adaptive resampling](@ref adaptiveresampling).
