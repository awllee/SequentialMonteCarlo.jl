# [Introduction](@id intro)

Sequential Monte Carlo (SMC) algorithms are defined in terms of an initial distribution $M_1$, a sequence of Markov transition kernels $M_{2}, \ldots, M_n$ and a sequence of non-negative potential functions $G_1, \ldots, G_n$.

These in turn define a sequence of unnormalized distributions, or measures, specified in [SMC Integrals](@ref smcintegrals). Integrals with respect to these measures are the principal quantities that an [SMC algorithm](@ref basicsmc) provides [approximations](@ref mainapprox) of.

The purpose of this package is to provide a very light [interface](@ref interface) to an SMC implementation with good multi-threaded performance. The main computationally intensive tasks in an SMC algorithm involve simulating from $M_1, \ldots, M_n$ and evaluating $G_1, \ldots, G_n$ a very large number of times: this package allows users to write these functions in Julia, which is a relatively high-level language, with little to no cost to performance. This is certainly an improvement over the very obsolete code written for Lee et al. (2010), which had good performance on GPUs but was almost impossible to make generally usable. It hopefully complements Murray (2013)'s [LibBi](http://libbi.org/) software, which is focused on Bayesian inference in the context of general state-space hidden Markov models, and which achieves genericity via the heavier interface of a custom modelling language and a high-performance back end.

## Which SMC algorithm?

The basic algorithm presented and implemented here was proposed by Stewart & McCarty Jr (1992), Gordon et al. (1993) and Kitagawa (1993). It is known as SMC with multinomial resampling and also often referred to as a Particle Filter. There are other SMC algorithms, some of which may be implemented in the future. The SMC algorithm with multinomial resampling is [theoretically very well understood](@ref maintheory), and has associated [variance estimators](@ref vepage) that are by-products of running the algorithm itself.

For readers interested in an SMC tutorial or methodological survey, two recent and complementary book chapters are Doucet & Johansen (2011) and Doucet & Lee (2018), which make some attempt at being comprehensive in various ways. The latter uses notation very similar to here, which is ultimately due to Pierre Del Moral. One minor notational difference here is the use of 1-indexing for sequences so as to be consistent with Julia's indexing. There are many other tutorials and surveys available.

Web pages maintained by [Arnaud Doucet](http://www.stats.ox.ac.uk/~doucet/smc_resources.html) and [Pierre del Moral](http://people.bordeaux.inria.fr/pierre.delmoral/simulinks.html) provide many important references to methodological, theoretical and applied work.

## Notation and assumptions

Let $(\mathsf{X}, \mathcal{X})$ be a measurable space. The initial distribution $M_1$ is a probability measure on this space, and the Markov kernels $M_{2}, \ldots, M_n$ evolve on this space. Every measure defined here, unless explicitly stated, will be a measure on this space. The non-negative functions $G_1, \ldots, G_n$ have domain $\mathsf{X}$. Every function $f : \mathsf{X} \rightarrow \mathbb{R}$ will be assumed to be measurable.

* *Weighting:* if $\mu$ is a measure and $g:\mathsf{X} \rightarrow \mathbb{R}_{+}$ a non-negative function, then $\mu \cdot g$ is the measure
```math
(\mu\cdot g)(A) := \int_{A}\mu({\rm d}x)g(x), \qquad A \in \mathcal{X}.
```
  One can think of $\mu \cdot g$ as $\mu$ **weighted** by $g$.

* *Mutation:* if $\mu$ is a measure and P is a Markov kernel then $\mu P$ is the measure defined by
```math
\mu P(A):= \int_{\mathsf{X}} \mu({\rm d}x) P(x,A), \qquad A \in \mathcal{X}.
```
  When $\mu$ is a probability measure, $\mu P(A) = \Pr(Y\in A)$ when $Y\sim P(X,\cdot)$ and $X \sim \mu$. One can think of $\mu P$ as $\mu$ **mutated** by the Markov transition kernel $P$.

* *Integration:* if $\mu$ is a measure and $f:\mathsf{X}\rightarrow\mathbb{R}$ then
```math
\mu(f) := \int_{\mathsf{X}}f(x)\mu({\rm d}x).
```

* We denote by $\delta_{x}$ the Dirac measure centred at x. Note that $\delta_{x}(f) = f(x)$.

* A random variable $X$ has a Categorical$(a_1,\ldots,a_{m})$ distribution if
```math
\Pr(X=i) = \frac{a_i}{\sum_{j=1}^m a_j} \mathbf{1}_{\{1,\ldots,m}\}(i).
```
