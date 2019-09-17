# Connection to hidden Markov models

Probably the most common use for SMC is in statistical applications involving hidden Markov models (HMMs), where the methodology originates. Indeed, there is a very simple correspondence between an HMM and an SMC model. The purpose of this part of the documentation is only to clarify how different SMC models can be associated with the same HMM.

An simple HMM can be described as a bivariate Markov chain $(X_1,Y_1), \ldots, (X_n, Y_n)$ as follows. $X_1, \ldots, X_n$ is a Markov chain determined by an initial distribution $\mu$ and Markov transition densities $f_2, \ldots, f_n$. For each $p \in \{1, \ldots, n\}$, $Y_p$ is conditionally independent of all other random variables given $X_p$ and $Y_p \mid (X_p = x)$ has density $g_p(x,  \cdot)$. Statistical inference for HMMs involve performing inference after observing $(y_1, \ldots, y_n)$ as a realization of $(Y_1, \ldots, Y_n)$.

A standard use of SMC for such an HMM is as follows. Let the observations $(y_1, \ldots, y_n)$ be given and choose $M_1 = \mu$, $M_p(x, dx') = f_p(x, x') dx'$ for $p \in \{2, \ldots, n\}$ and $G_p(x) = g_p(x, y_p)$ for $p \in \{1, \ldots, n\}$. Notice that the potential function $G_p$ is defined using the observation $y_p$. Then it can be deduced that $\hat{Z}_p$ is the *marginal likelihood* associated with the first $p$ observations, $\eta_p$ is the *predictive* distribution of $X_p$ given the first $p-1$ observations and $\hat{\eta}_p$ is the *filtering* distribution of $X_p$ given the first $p$ observations. That is,
```math
X_p \mid \\{(Y_0, \ldots, Y_{p-1}) = (y_0, \ldots, y_{p-1})\\} \sim \eta_p
```
and
```math
X_p \mid \\{(Y_0, \ldots, Y_{p}) = (y_0, \ldots, y_{p})\\} \sim \hat{\eta}_p.
```

However, this is only one way to map the HMM into an SMC Model. We extend the HMM's state space and define $M_1(d(x_{\rm prev}, x)) = \mu(dx) \delta_x(dx_{\rm prev})$, and
```math
M_p((x_{\rm prev},x), d(x_{\rm prev}',x')) = q_p(x, x') \delta_x(dx_{\rm prev}'), \qquad p \in \{2, \ldots, n\},
```
where $q_p$ is such that $q_p(x, x') > 0$ whenever $f_p(x, x') > 0$. Then one can define
```math
G_p(x_{\rm prev}, x) = \frac{f_p(x_{\rm prev}, x)}{q_p(x_{\rm prev}, x)} g(x, y_p), \qquad p \in \{1, \ldots, n\}
```
and it is still the case that $\hat{Z}_p$ is the marginal likelihood associated with the first $p$ observations, $\eta_p$ is the predictive distribution of $X_p$ given the first $p-1$ observations and $\hat{\eta}_p$ is the filtering distribution of $X_p$ given the first $p$ observations. An example of using different SMC models for a given HMM is provided in [SMCExamples: comparison of Linear Gaussian models](https://github.com/awllee/SMCExamples.jl/blob/master/demo/lgComparisonDemo.jl) ; the model associated with the "locally optimal proposal" is defined [here](https://github.com/awllee/SMCExamples.jl/blob/master/src/lgLOPModel.jl).
