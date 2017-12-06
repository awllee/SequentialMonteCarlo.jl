var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Contents",
    "title": "Contents",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SequentialMonteCarlo.jl-1",
    "page": "Contents",
    "title": "SequentialMonteCarlo.jl",
    "category": "section",
    "text": "A light interface to serial and multi-threaded Sequential Monte CarloPages = [ \"intro.md\", \"smcintegrals.md\", \"smcalgorithm.md\", \"smctheory.md\", \"smcve.md\", \"smcadaptive.md\", \"csmc.md\", \"impl.md\", \"smcinterface.md\", \"guide.md\", \"refs.md\"]"
},

{
    "location": "intro.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "intro.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "Sequential Monte Carlo (SMC) algorithms are defined in terms of an initial distribution M_1, a sequence of Markov transition kernels M_2 ldots M_n and a sequence of non-negative potential functions G_1 ldots G_n.These in turn define a sequence of unnormalized distributions, or measures, specified in SMC Integrals. Integrals with respect to these measures are the principal quantities that an SMC algorithm provides approximations of.The purpose of this package is to provide a very light interface to an SMC implementation with good multi-threaded performance. The main computationally intensive tasks in an SMC algorithm involve simulating from M_1 ldots M_n and evaluating G_1 ldots G_n a very large number of times: this package allows users to write these functions in Julia, which is a relatively high-level language, with little to no cost to performance. This is certainly an improvement over the very obsolete code written for Lee et al. (2010), which had good performance on GPUs but was almost impossible to make generally usable. It hopefully complements Murray (2013)'s LibBi software, which is focused on Bayesian inference in the context of general state-space hidden Markov models, and which achieves genericity via the heavier interface of a custom modelling language and a high-performance back end."
},

{
    "location": "intro.html#Which-SMC-algorithm?-1",
    "page": "Introduction",
    "title": "Which SMC algorithm?",
    "category": "section",
    "text": "The basic algorithm presented and implemented here was proposed by Stewart & McCarty Jr (1992), Gordon et al. (1993) and Kitagawa (1993). It is known as SMC with multinomial resampling and also often referred to as a Particle Filter. There are other SMC algorithms, some of which may be implemented in the future. The SMC algorithm with multinomial resampling is theoretically very well understood, and has associated variance estimators that are by-products of running the algorithm itself.For readers interested in an SMC tutorial or methodological survey, two recent and complementary book chapters are Doucet & Johansen (2011) and Doucet & Lee (2018), which make some attempt at being comprehensive in various ways. The latter uses notation very similar to here, which is ultimately due to Pierre Del Moral. One minor notational difference here is the use of 1-indexing for sequences so as to be consistent with Julia's indexing. There are many other tutorials and surveys available.Web pages maintained by Arnaud Doucet and Pierre del Moral provide many important references to methodological, theoretical and applied work."
},

{
    "location": "intro.html#Notation-and-assumptions-1",
    "page": "Introduction",
    "title": "Notation and assumptions",
    "category": "section",
    "text": "Let (mathsfX mathcalX) be a measurable space. The initial distribution M_1 is a probability measure on this space, and the Markov kernels M_2 ldots M_n evolve on this space. Every measure defined here, unless explicitly stated, will be a measure on this space. The non-negative functions G_1 ldots G_n have domain mathsfX. Every function f  mathsfX rightarrow mathbbR will be assumed to be measurable.Weighting: if mu is a measure and gmathsfX rightarrow mathbbR_+ a non-negative function, then mu cdot g is the measure \\[ (\\mu\\cdot g)(A) := \\int_{A}\\mu({\\rm d}x)g(x), \\qquad A \\in \\mathcal{X}. \\] One can think of mu cdot g as mu weighted by g.\nMutation: if mu is a measure and P is a Markov kernel then mu P is the measure defined by \\[ \\mu P(A):= \\int_{\\mathsf{X}} \\mu({\\rm d}x) P(x,A), \\qquad A \\in \\mathcal{X}. \\] When mu is a probability measure, mu P(A) = Pr(Yin A) when Ysim P(Xcdot) and X sim mu. One can think of mu P as mu mutated by the Markov transition kernel P.\nIntegration: if mu is a measure and fmathsfXrightarrowmathbbR then \\[ \\mu(f) := \\int_{\\mathsf{X}}f(x)\\mu({\\rm d}x). \\]\nWe denote by delta_x the Dirac measure centred at x. Note that delta_x(f) = f(x).\nA random variable X has a Categorical(a_1ldotsa_m) distribution if \\[ \\Pr(X=i) = \\frac{a_i}{\\sum_{j=1}^m a_j} \\mathbf{1}_{\\{1,\\ldots,m\\}}(i). \\]"
},

{
    "location": "smcintegrals.html#",
    "page": "SMC integrals",
    "title": "SMC integrals",
    "category": "page",
    "text": ""
},

{
    "location": "smcintegrals.html#smcintegrals-1",
    "page": "SMC integrals",
    "title": "Integrals approximated by SMC",
    "category": "section",
    "text": ""
},

{
    "location": "smcintegrals.html#The-measures-and-integrals-1",
    "page": "SMC integrals",
    "title": "The measures and integrals",
    "category": "section",
    "text": "We define a sequence of measures gamma_1 ldots gamma_n by gamma_1 = M_1 and for p in 2ldotsn, \\[ \\gamma_p(A) := \\int_{\\mathsf{X}^p} \\mathbf{1}_A(x_p) M_1({\\rm d}x_1) \\prod_{q=2}^p G_{q-1}(x_{q-1}) M_q(x_{q-1}, {\\rm d}x_q), \\qquad A \\in \\mathcal{X}. \\] We also define a second sequence of measures hatgamma_1ldotshatgamma_n by hatgamma_p = gamma_p cdot G_p for p in 1ldotsn.Letting 1 denote the constant function x mapsto 1, the normalizing constants associated with each gamma_p and hatgamma_p are gamma_p(1) and hatgamma_p(1), respectively. We define \\[ Z_p := \\gamma_p(1), \\qquad \\hat{Z}_p := \\hat{\\gamma}_p(1), \\qquad p\\in \\{1,\\ldots,n\\}. \\] This allows us to define their normalized, probability measure counterparts \\[ \\eta_p := \\gamma_p / Z_p, \\qquad \\hat{\\eta}_p := \\hat{\\gamma}_p / \\hat{Z}_p, \\qquad p \\in \\{1,\\ldots,n\\}. \\]We are now able to state which integrals SMC algorithms are built to approximate: these are simply gamma_p(f), hatgamma_p(f), eta_p(f) or hateta_p(f) for some p in 1ldotsn and some function fmathsfXrightarrow mathbbR.We note from their definitions that Z_1 = 1 and Z_p = hatZ_p-1 for p in 2ldotsn. It follows that any of the above measures can be written in terms of the sequences hatZ_1 ldots hatZ_n, eta_1 ldots eta_n, and hateta_1 ldots hateta_n. For example, gamma_p = hatZ_p-1 eta_p for p in 2 ldots n. Therefore, it is sufficient to define SMC approximations of hatZ_p, eta_p and hateta_p for p in 1 ldots n ."
},

{
    "location": "smcintegrals.html#recursivedefmeasures-1",
    "page": "SMC integrals",
    "title": "Recursive definition",
    "category": "section",
    "text": "We can alternatively define the sequence of measures mentioned above recursively by eta_1=gamma_1=M_1, \\[ \\eta_p \\propto \\gamma_p := (\\gamma_{p-1} \\cdot G_{p-1}) M_p, \\qquad p \\in \\{2, \\ldots, n\\}, \\] and for pin 1ldotsn, hateta_p propto hatgamma_p = gamma_p cdot G_p. It follows that one can write hatZ_1 = eta_1(G_1) and hatZ_p = hatZ_p-1 eta_p(G_p) for p in 2 ldots n.We observe that gamma_p=hatgamma_p-1M_p, and this suggests the following diagram \\[ M_1 = \\gamma_1 \\overset{G_1}{\\longrightarrow} \\hat{\\gamma}_1 \\overset{M_{2}}{\\longrightarrow} \\gamma_{2} \\longrightarrow \\cdots \\longrightarrow \\hat{\\gamma}_{p-1} \\overset{M_p}{\\longrightarrow} \\gamma_p \\overset{G_p}{\\longrightarrow} \\hat{\\gamma}_p \\longrightarrow \\cdots \\longrightarrow \\gamma_n \\overset{G_n}{\\longrightarrow} \\hat{\\gamma}_n, \\] where measures are obtained by weighting or mutation. Alternatively, we can view the sequence of probability measures using the same diagram, i.e. \\[ M_1 = \\eta_1 \\overset{G_1}{\\longrightarrow} \\hat{\\eta}_1 \\overset{M_{2}}{\\longrightarrow} \\eta_{2} \\longrightarrow \\cdots \\longrightarrow \\hat{\\eta}_{p-1} \\overset{M_p}{\\longrightarrow} \\eta_p \\overset{G_p}{\\longrightarrow} \\hat{\\eta}_p \\longrightarrow \\cdots \\longrightarrow \\eta_n \\overset{G_n}{\\longrightarrow} \\hat{\\eta}_n, \\] but where the weighting steps also involve renormalization: eta_p oversetG_plongrightarrow hateta_p means hateta_p = (eta_p cdot G_p)  eta_p(G_p).The measures hatgamma_p and hateta_p are often referred to as the \"updated\" versions of gamma_p and eta_p, respectively."
},

{
    "location": "smcalgorithm.html#",
    "page": "SMC algorithm",
    "title": "SMC algorithm",
    "category": "page",
    "text": ""
},

{
    "location": "smcalgorithm.html#SMC-algorithm-and-particle-approximations-1",
    "page": "SMC algorithm",
    "title": "SMC algorithm and particle approximations",
    "category": "section",
    "text": ""
},

{
    "location": "smcalgorithm.html#basicsmc-1",
    "page": "SMC algorithm",
    "title": "SMC algorithm",
    "category": "section",
    "text": "The implemented algorithm is as follows. The random variables zeta_p^1 ldots zeta_p^N are typically referred to as particles, and the approximations they define are referred to as particle approximations.Sample zeta_1^i oversetmathrmiidsim M_1 for i in 1ldotsN.\nFor p=2ldotsn:\nsample a vector A_p-1^1ldotsA_p-1^N of i.i.d. rm Categorical(G_p-1(zeta_p-1^1)ldotsG_p-1(zeta_p-1^N)) random variables in increasing order.\nsample zeta_p^i oversetmathrmind sim M_p(zeta_p-1^A_p-1^i cdot) for i in 1ldotsN.This algorithm is referred to as SMC with multinomial resampling because the vector of numbers of offspring of the time p-1 particles follows a multinomial distribution."
},

{
    "location": "smcalgorithm.html#mainapprox-1",
    "page": "SMC algorithm",
    "title": "Main particle approximations",
    "category": "section",
    "text": "For some p in 1ldotsn, the particle approximation of eta_p(f) is obtained by constructing a particle approximation eta_p^N of eta_p and then calculating eta_p^N(f). Specifically, we define \\[ \\eta_p^N := \\frac{1}{N} \\sum_{i=1}^N \\delta_{\\zeta_p^i}, \\qquad p \\in \\{1,\\ldots,n\\} \\] and the particle approximation of eta_p(f) is then simply \\[ \\eta_p^N(f) = \\int f(x) \\eta_p^N({\\rm d}x) = \\frac{1}{N} \\sum_{i=1}^N f(\\zeta_p^i). \\]Similarly, we define \\[ \\hat{\\eta}_p^N:= \\frac{\\eta_p^N \\cdot G_p}{\\eta_p^N(G_p)}, \\qquad p \\in \\{1, \\ldots, n\\} \\] so that the particle approximation of hateta_p(f) is \\[ \\hat{\\eta}^N_p(f) = \\int f(x) \\hat{\\eta}_p^N({\\rm d}x) \\frac{\\sum_{i=1}^N G_p(\\zeta_p^i) f(\\zeta_p^i)}{\\sum_{i=1}^N G_p(\\zeta_p^i)}. \\]Finally, the particle approximations of hatZ_1ldotshatZ_n are defined by hatZ_1^N=eta_1^N(G_1) and for p in 2ldotsn, hatZ_p^N=hatZ_p-1^Neta_p^N(G_p)."
},

{
    "location": "smcalgorithm.html#Genealogical-structure-1",
    "page": "SMC algorithm",
    "title": "Genealogical structure",
    "category": "section",
    "text": "The particle system defined by the SMC algorithm has a genealogical structure, induced by the ancestor indices. We define Eve indices such that E_p^i is the index of the time 1 ancestor of zeta_p^i. That is, the Eve indices satisfy E_1^i = i for i in 1 ldots N and \\[ E_p^i = E_{p-1}^{A_{p-1}^i}, \\qquad p \\in \\{2, \\ldots, n\\}, \\quad i \\in \\{1, \\ldots, N\\}. \\] These can be used to quantify and visualize the path degeneracy phenomenon associated with the SMC algorithm. Since the ancestor indices are in sorted order, so too are the Eve indices at each time."
},

{
    "location": "smcalgorithm.html#Recursive-definition-1",
    "page": "SMC algorithm",
    "title": "Recursive definition",
    "category": "section",
    "text": "The SMC algorithm is intimately connected to the recursive definition of the SMC measures themselves.Let us define the sampling operator S^N to be the random map that takes a measure mu and outputs a random measure \\[ \\mu^N = S^N \\mu = \\frac{1}{N} \\sum_{i=1}^N \\delta_{X^i}, \\qquad \\text{where } X^i \\overset{\\mathrm{i.i.d.}}{\\sim} \\mu. \\]Then we have \\[ \\eta^N_1 = S^N M_1, \\qquad \\hat{\\eta}_1^N = \\frac {\\eta^N_1 \\cdot G_1} {\\eta^N_1(G_1)}, \\] and for p = 2 ldots n, \\[ \\eta^N_p = S^N (\\hat{\\eta}_{p-1}^N M_p), \\qquad \\hat{\\eta}_p^N = \\frac {\\eta^N_p \\cdot G_p} {\\eta^N_p(G_p)}. \\]This mirrors almost exactly the recursive definitions for the SMC measures, where one has instead \\[ \\eta_1 = M_1, \\qquad \\hat{\\eta}_1 = \\frac {\\eta_1 \\cdot G_1} {\\eta_1(G_1)}, \\] and for p = 2 ldots n, \\[ \\eta_p = \\hat{\\eta}_{p-1} M_p, \\qquad \\hat{\\eta}_p = \\frac {\\eta_p \\cdot G_p} {\\eta_p(G_p)}. \\]We notice also that hatZ^N_p = hatZ^N_p-1 eta^N_p(G_p) while hatZ_p = hatZ_p-1 eta_p(G_p) for p in 2 ldots n."
},

{
    "location": "smctheory.html#",
    "page": "Theoretical properties",
    "title": "Theoretical properties",
    "category": "page",
    "text": ""
},

{
    "location": "smctheory.html#maintheory-1",
    "page": "Theoretical properties",
    "title": "Theoretical properties of particle approximations",
    "category": "section",
    "text": "All convergence results below involve taking limits as Nrightarrowinfty. Convergence almost surely is denoted oversetasrightarrow, convergence in probability is denoted oversetPrightarrow and weak convergence or convergence in distribution/law is denoted oversetLrightarrow.The following Theorem provides some justification for the use of the main particle approximations; these now classical results can be deduced from various results of Del Moral (2004). We present the results only for the \"hatted\" or \"updated\" quantities to avoid repetition; analogous hold for their \"unhatted\" counterparts with a different sequence of maps sigma_1^2 ldots sigma_n^2.Theorem [Del Moral, 2004]. Let the potential functions G_1 ldots G_n be bounded and strictly positive. There exist maps hatsigma_1^2 ldots hatsigma_n^2 such that the following hold for an arbitrary, bounded f:Lack-of-bias: mathbbE left hatgamma_p^N(f) right = hatgamma_p(f) for all N geq 1.\nConsistency: hatgamma_p^N(f)oversetasrightarrowhatgamma_p(f) and hateta_p^N(f)oversetasrightarrowhateta_p(f).\nAsymptotic variance and mean-squared error (MSE): \\[ N {\\rm var} \\left\\{ \\hat{\\gamma}_p^N(f) / \\hat{\\gamma}_p(1) \\right\\} \\rightarrow \\hat{\\sigma}_p^2(f), \\] and \\[ N \\mathbb{E} \\left[ \\left\\{ \\hat{\\eta}_p^N(f) - \\hat{\\eta}_p(f) \\right\\} ^2 \\right] \\rightarrow \\hat{\\sigma}_p^2(f-\\hat{\\eta}_p(f)). \\]\nCentral Limit Theorems: \\[ \\sqrt{N} \\left( \\hat{\\gamma}_p^N(f) / \\hat{\\gamma}_p(1) - \\eta_p^N(f) \\right) \\overset{L}{\\rightarrow} N(0,\\hat{\\sigma}_p^2(f)), \\] and \\[ \\sqrt{N} \\left( \\hat{\\eta}_p^N(f) - \\eta_p^N(f) \\right) \\overset{L}{\\rightarrow} N(0,\\hat{\\sigma}_p^2(f - \\hat{\\eta}_p(f))). \\]"
},

{
    "location": "smctheory.html#Note-on-the-sorted-ancestor-indices-1",
    "page": "Theoretical properties",
    "title": "Note on the sorted ancestor indices",
    "category": "section",
    "text": "The theoretical results above are typically proven for an algorithm that differs very slightly from the SMC algorithm implemented here. In particular, one would usually analyze the algorithm by considering A_p-1^1ldotsA_p-1^N to be i.i.d. rm Categorical (G_p-1(zeta_p-1^1) ldots G_p-1(zeta_p-1^N)) random variables rather than being in sorted order.The laws of the approximations hatZ_p^N, eta_p^N(f) and hateta_p^N(f), however, are invariant to permutations of the indices of the ancestors A_p-1^1ldotsA_p-1^N in the algorithm. Therefore, the sorting of the ancestor indices may be regarded as an implementation issue that does not affect the particle approximations themselves."
},

{
    "location": "smcve.html#",
    "page": "Variance estimators",
    "title": "Variance estimators",
    "category": "page",
    "text": ""
},

{
    "location": "smcve.html#vepage-1",
    "page": "Variance estimators",
    "title": "Variance estimators",
    "category": "section",
    "text": ""
},

{
    "location": "smcve.html#Asymptotic-variances-1",
    "page": "Variance estimators",
    "title": "Asymptotic variances",
    "category": "section",
    "text": "We recall that the particle approximations have asymptotic variances associated with particular test functions f. Letting n be arbitrary to avoid additional indexing notation, the values sigma_n^2(f) and hatsigma_n^2(f) are typically not possible to calculate exactly. They can be decomposed into sums as follows: \\[ \\sigma_n^2(f) = \\sum_{p=1}^n v_{p,n}(f), \\qquad \\hat{\\sigma}_n^2(f) = \\sum_{p=1}^n \\hat{v}_{p,n}(f).\\] The quantities v_pn(f) and hatv_pn(f) also cannot be computed exactly in general, but this decomposition can in some circumstances shed some light on the nature of the approximation errors."
},

{
    "location": "smcve.html#The-estimators-1",
    "page": "Variance estimators",
    "title": "The estimators",
    "category": "section",
    "text": "Most Monte Carlo approximations are accompanied by Monte Carlo approximations of their asymptotic variance or MSE. Lee and Whiteley (2015) defined V_p^N(f), computable using the simulated particle system, allowing one to approximate asymptotic variances of interest.approximation approximation of\nV_p^N(f) rm var left  gamma_p^N(f)gamma_p(1) right \nV_p^N(f-eta_p^N(f)) mathbbE left left eta_p^N(f) - eta_p(f) right^2 right\nhatV_p^N(f) rm var left hatgamma_p^N(f)  hatgamma_p(1) right\nhatV_p^N(f-hateta_p^N(f)) mathbbE left left hateta_p^N(f) - hateta_p(f) right ^2rightThese approximations have the following theoretical justification, particularly for their approximation of the asymptotic variance and MSE maps defined earlier. Again, we present the results only for the \"updated\" quantities to avoid repetition; analogous results hold for their \"unhatted\" counterparts."
},

{
    "location": "smcve.html#Theoretical-justification-1",
    "page": "Variance estimators",
    "title": "Theoretical justification",
    "category": "section",
    "text": "The variance approximations have consistency properties, and some of them have lack-of-bias properties.Theorem [Lee & Whiteley, 2015]. Let the potential functions G_1 ldots G_n be bounded and strictly positive. The following hold for an arbitrary, bounded f:Lack-of-bias: mathbbEleft(hatZ_p^N)^2hatV_p^N(f)right=rm varleft hatgamma_p(f)right for all N geq 1.\nConsistency: NhatV_p^N(f) oversetPrightarrow hatsigma_p^2(f) and N hatV_p^N(f-hateta_p^N(f)) oversetPrightarrow hatsigma_p^2(f - hateta_p(f)).Remark. A consistent estimator of the asymptotic MSE of eta_p(f) or hateta_p(f) was proposed in Chan & Lai (2013). This is very similar to that of Lee & Whiteley (2015) but does not satisfy the lack-of-bias property.Lee & Whiteley (2015) also define approximations v_pn^N(f) and hatv_pn^N(f) of v_pn(f) and hatv_pn(f), respectively. Their sums v_n^N(f)=sum_p=1^nv_pn^N(f) and hatv_n^N(f)=sum_p=1^nhatv_pn^N(f) can also be used as alternative approximations of sigma_n^2(f) and hatsigma_n^2(f)Theorem [Lee & Whiteley, 2015]. Let the potential functions G_1 ldots G_n be bounded and strictly positive. The following hold for an arbitrary, bounded f:Lack-of-bias: mathbbEleftleft(hatZ_n^Nright)^2hatv_pn^N(f)right=hatZ_n^2hatsigma_n^2(f) for all Ngeq1.\nConsistency: hatv_pn^N(f)oversetPrightarrowhatv_pn(f) and hatv_pn^N(f)(f-hateta_n^N(f))oversetPrightarrowhatv_pn(f-hateta_n(f)).\nLack-of-bias: mathbbEleftleft(hatZ_n^Nright)^2hatv_n^N(f)right=hatZ_n^2hatsigma_n^2(f) for all Ngeq1.\nConsistency: hatv_n^N(f)oversetPrightarrowhatsigma_n^2(f) and hatv_n^N(f-hateta_n^N(f))oversetPrightarrowhatsigma_n^2(f-hateta_n(f))."
},

{
    "location": "smcadaptive.html#",
    "page": "Adaptive resampling",
    "title": "Adaptive resampling",
    "category": "page",
    "text": ""
},

{
    "location": "smcadaptive.html#adaptiveresampling-1",
    "page": "Adaptive resampling",
    "title": "SMC with adaptive resampling",
    "category": "section",
    "text": "There is a tuning parameter associated with the SequentialMonteCarlo.smc algorithm that results in an adaptive version of the SMC Algorithm. This is the essThreshold parameter, whose use was proposed by Kong et al. (1994) and Liu & Chen (1995). We represent this parameter as tau below."
},

{
    "location": "smcadaptive.html#(Relative)-effective-sample-size-1",
    "page": "Adaptive resampling",
    "title": "(Relative) effective sample size",
    "category": "section",
    "text": "The adaptive SMC algorithm essentially involves choosing (A_p-1^1 ldots A_p-1^N) = left(1 ldots N right) when the weights associated with particles zeta_p-1^1 ldots zeta_p-1^N have a relative effective sample size exceeding tau. This is sometimes referred to as \"not resampling\". The relative effective sample size of a collection of weights is defined as \\[ {\\rm rESS}(w_{1}, \\ldots, w_{N}) := \\frac{\\left(\\frac{1}{N}\\sum_{i=1}^N w_{i} \\right)^{2}}{\\frac{1}{N} \\sum_{i=1}^N w_{i}^{2}}. \\] This function is invariant to rescaling of all of the weights by a common constant.In the SMC Algorithm we can view the weights associated with the particles zeta_p-1^1 ldots zeta_p-1^N as being W_p-1^1 ldots W_p-1^N where W_p-1^i propto G_p-1(zeta_p-1^i), and their relative effective sample size is rm rESS(W_p-1^1 ldots W_p-1^N).In the SMC with adaptive resampling algorithm, the weight of a particle may be proportional to a product of many potential function values, depending on when the most recent normal resampling step was. In particular, whenever the standard resampling step is not taken, particles first inherit the weights of their ancestors and then multiply them by their own potential function values. The algorithm is as follows:"
},

{
    "location": "smcadaptive.html#Adaptive-resampling-algorithm-1",
    "page": "Adaptive resampling",
    "title": "Adaptive resampling algorithm",
    "category": "section",
    "text": "Sample zeta_1^i oversetmathrmiidsim M_1 and compute W_1^i propto G_1(zeta_1^i) for iin 1 ldots N.\nFor p=2ldotsn:\ncompute mathcalE_p-1^N = rm rESS(W_p-1^1ldotsW_p-1^N).\nif p = n or mathcalE_p-1^N leq tau set R_p-1 leftarrow 1; otherwise set R_p-1 leftarrow 0.\nif R_p-1 = 1, sample a vector A_p-1^1 ldots A_p-1^N of i.i.d. rm Categorical(W_p-1^1 ldots W_p-1^N) random variables in increasing order; otherwise set (A_p-1^1 ldots A_p-1^N) = left(1 ldots N right).\nsample zeta_p^i oversetmathrmindsim M_p(zeta_p-1^A_p-1^i cdot) and compute W_p^i propto left(W_p-1^i right)^mathbbI(R_p-1 = 0) G_p(zeta_p^i) for iin 1ldotsN.Note that the time n-1 particles are always resampled, so that eta_n^N is always an unweighted particle approximation of eta_n. This is primarily an implementation detail."
},

{
    "location": "smcadaptive.html#Particle-approximations-1",
    "page": "Adaptive resampling",
    "title": "Particle approximations",
    "category": "section",
    "text": "We define \\[ \\eta_{p}^N \\propto \\sum_{i=1}^N \\left(W_{p-1}^{i} \\right)^{ \\mathbb{I}(R_{p-1} = 0)} \\delta_{\\zeta_{p}^{i}}, \\qquad p \\in \\{1,\\ldots,n\\}, \\] and \\[ \\hat{\\eta}_{p}^N \\propto \\sum_{i=1}^N W_{p}^{i} \\delta_{\\zeta_{p}^{i}}, \\qquad p \\in \\{1,\\ldots,n\\}. \\]Appropriate approximations of hatZ_1^N ldots hatZ_n^N are also well-defined, but tedious to display.These particle approximations all enjoy the same theoretical properties stated for the standard SMC Algorithm, although the asymptotic variance maps sigma^2_p are generally different. In fact, adaptive resampling improves in certain scenarios the quality of SMC approximations; its effects have been studied theoretically by Del Moral et al. (2010) and Whiteley et al. (2016)."
},

{
    "location": "smcadaptive.html#Variance-estimation-1",
    "page": "Adaptive resampling",
    "title": "Variance estimation",
    "category": "section",
    "text": "All of the methods described in Variance estimators can be run on output from the SMC with adaptive resampling algorithm.One should be aware that the length of the vector returned by SequentialMonteCarlo.vpns will be of length m = 1 + sum_i=1^n-1 R_i leq n. This is a consequence of the fact that resampling only at certain times can be viewed as running an SMC algorithm with modified Markov kernels and potential functions defined on a more complicated state space –- the details are not presented here –- so that the number of terms in the asymptotic variance decomposition is m and not (necessarily) n."
},

{
    "location": "csmc.html#",
    "page": "Conditional SMC",
    "title": "Conditional SMC",
    "category": "page",
    "text": ""
},

{
    "location": "csmc.html#Conditional-SMC-1",
    "page": "Conditional SMC",
    "title": "Conditional SMC",
    "category": "section",
    "text": "One relatively recent subtle modification of the SMC algorithm was proposed by Andrieu et al. (2010), and known as conditional SMC (cSMC). The algorithm is as follows:Input: reference path (x_1^rm refldotsx_n^rm ref)Set zeta_1^1=x_1^rm ref and sample zeta_1^ioversetmathrmiidsimM_1 for i in 2 ldots N .\nFor p=2ldotsn:\nset A_p-1^1=1 and zeta_p^1=x_p^rm ref.\nsample a vector A_p-1^2ldotsA_p-1^N of i.i.d. rm Categorical(G_p-1(zeta_p-1^1)ldotsG_p-1(zeta_p-1^N)) random variables in increasing order.\nsample zeta_p^i oversetmathrmindsim M_p(zeta_p-1^A_p-1^i cdot) for i in 2 ldots N .\nSample K_n sim rm Categorical(G_n(zeta_n^1)ldotsG_n(zeta_n^N)) and for p=n-1ldots1 set K_p=A_p^K_p+1.\nOutput the path (zeta_1^K_1ldotszeta_n^K_n).cSMC is of special interest as it defines a Markov kernel that is ergodic, for N geq 2, with invariant probability measure given by \\[ \\hat{\\pi}(A)=\\hat{Z}_{n}^{-1}\\int_{A}G_{n}(x_{n})M_{1}({\\rm d}x_{1})\\prod_{p=2}^{n}G_{p-1}(x_{p-1})M_{p}(x_{p-1},{\\rm d}x_{p}), \\] for A in mathcalX^otimesn.That is, the Markov kernel P_N(x^rm ref cdot) defined by calling the cSMC algorithm has invariant probability measure hatpi. Theoretical results, in particular concerning the effect of N on the rate of convergence of P_N, have been developed by Chopin & Singh (2015), Andrieu et al. (2018) and Lindsten et al. (2015).Conditional SMC may also be used with adaptive resampling as detailed in Adaptive resampling."
},

{
    "location": "impl.html#",
    "page": "Implementation notes",
    "title": "Implementation notes",
    "category": "page",
    "text": ""
},

{
    "location": "impl.html#Implementation-notes-1",
    "page": "Implementation notes",
    "title": "Implementation notes",
    "category": "section",
    "text": "The multi-threaded implementation of SMC is fairly straightforward. The attempt to efficiently utilize multiple cores is largely focused around low-overhead resampling / selection, since the remaining steps are embarrassingly parallel.There is also some attention paid to Julia's threading interface and, for both parallel and serial execution, the need to avoid dynamic memory allocations. All code for the main package was written using Julia's core and base functionality.Generating the ancestor indices in sorted order is accomplished primarily through use of the uniform spacings method for generating a sequence of sorted uniform random variates developed by Lurie & Hartley (1972) and described by Devroye (1986, p. 214). In multi-threaded code, the multinomial variate giving the number of ancestor indices each thread should simulate using a thread-specific local vector of particle weights is simulated by using a combination of inversion sampling and a a Julia implementation of the btrd algorithm of Hörmann (1993), which simulates binomially distributed variates in expected mathcalO(1) time. Simulating a multinomial random variate is accomplished by simulating a sequence of appropriate Binomial random variates. The code-generating function for copying particles was adapted from a suggestion on Julia Discourse by Greg Plowman.For the demonstration code, use of the StaticArrays.jl package provides dramatic improvements for multivariate particles, and is highly recommended for fixed-size vector components of particles. It is possible to use the counter-based RNGs in the RandomNumbers.jl package in place of the default random number generator; this was not working on the Julia-0.7 master branch at the time of development, and is slightly slower than the MersenneTwister RNG provided by Julia's base."
},

{
    "location": "smcinterface.html#",
    "page": "SMC interface",
    "title": "SMC interface",
    "category": "page",
    "text": ""
},

{
    "location": "smcinterface.html#interface-1",
    "page": "SMC interface",
    "title": "Interface",
    "category": "section",
    "text": ""
},

{
    "location": "smcinterface.html#Specifying-an-SMC-model-1",
    "page": "SMC interface",
    "title": "Specifying an SMC model",
    "category": "section",
    "text": "We recall from the Introduction that the SMC algorithm is defined in terms of M_1 ldots M_n and G_1 ldots G_n.A model of type SMCModel can be created by callingmodel = SMCModel(M!, lG, maxn::Int64, Particle, ParticleScratch)where Particle and ParticleScratch are user-defined types, M! is a void function with argumentsM!(newParticle::Particle, rng::SMCRNG, p::Int64, particle::Particle,\n  scratch::ParticleScratch)and lG is a function returning a Float64 and has argumentslG(p::Int64, particle::Particle, scratch::ParticleScratch)There is a correspondence between the function M! and M_1 ldots M_n: calling M!(x', rng, p, x, scratch), should make x a realization of a sample from M_p(x cdot) with the convention that M_1(xcdot) = M_1(cdot) for any x. Similarly, lG and G_1 ldots G_n correspond in that lG(p,x) $ = \\log G_p(x)$. Logarithms are used to avoid numerical issues. maxn is the maximum value of n for which M! and lG are well-defined; users may choose to run the SMC algorithm for any integer value of n up to and including maxn.The types Particle and ParticleScratch must have constructors that take no arguments. One may choose ParticleScratch = Void, in which case nothing will be passed to M! and lG. Using scratch space is optional but can significantly improve performance in certain scenarios; it provides a mechanism for users to avoid dynamic memory allocations in M! and/or lG. This scratch space will be used by every particle associated with a given thread. A thread-specific pseudo-random number generator (RNG) rng will be passed to the M! function by the algorithm, and should be used in lieu of Julia's global RNG."
},

{
    "location": "smcinterface.html#Running-the-SMC-algorithm-1",
    "page": "SMC interface",
    "title": "Running the SMC algorithm",
    "category": "section",
    "text": "The SMC algorithm is run by callingsmc!(model::SMCModel, smcio::SMCIO)The second argument, smcio, is a struct containing inputs and outputs for the smc! algorithm. It is straightforward to construct and involves specifying the number of particles N, the number of iterations n, the number of threads nthreads, whether the entire history of the particle system should be recorded fullOutput, and an effective sample size threshold essThreshold that is explained on the adaptive resampling page, and can be safely ignored on a first reading. smcio can be created by callingsmcio = SMCIO{model.particle, model.pScratch}(N::Int64, n::Int64,\n  nthreads::Int64, fullOutput::Bool, essThreshold::Float64 = 2.0)Currently N must be an integer multiple of nthreads; this may change in the future."
},

{
    "location": "smcinterface.html#Extracting-the-main-outputs-1",
    "page": "SMC interface",
    "title": "Extracting the main outputs",
    "category": "section",
    "text": "A vector of approximations (loghatZ_1^N ldots loghatZ_n^N) is stored in smcio.logZhats, which can be used in conjunction with SequentialMonteCarlo.eta to produce the approximations gamma_p^N(f) or hatgamma_p^N(f).One can extract the approximation eta_p^N(f) or hateta_p^N(f) by callingSequentialMonteCarlo.eta(smcio, f, hat::Bool, p::Int64)with hat determining which approximation is returned. Calling this function with p < smcio.n requires smcio.fullOutput = true.One can extract the approximations (eta^N_p(f) geq 0 log gamma^N_p(f)) or (hateta^N_p(f) geq 0 log hatgamma^N_p(f)) by callingslgamma(smcio, f, hat::Bool, p::Int64)with hat determining which approximation is returned. Calling this function with p < smcio.n requires smcio.fullOutput = true."
},

{
    "location": "smcinterface.html#Variance-estimators-1",
    "page": "SMC interface",
    "title": "Variance estimators",
    "category": "section",
    "text": "The following functions relate to the approximations detailed in the variance estimators page.The functionSequentialMonteCarlo.V(smcio, f, hat::Bool, centred::Bool, p::Int64)returns:hat centred approximation approximation of\nfalse false V_p^N(f) rm var left  gamma_p^N(f)gamma_p(1) right \nfalse true V_p^N(f-eta_p^N(f)) mathbbEleftleft eta_p^N(f)-eta_p(f)right ^2right\ntrue false hatV_p^N(f) rm varleft hatgamma_p^N(f)hatgamma_p(1)right\ntrue true hatV_p^N(f-hateta_p^N(f)) mathbbEleftleft hateta_p^N(f)-hateta_p(f)right ^2rightCalling SequentialMonteCarlo.V with p < smcio.n requires smcio.fullOutput = true.A vector of the quantities hatV_1^N(1)ldotshatV_n^N(1) is stored in smcio.Vhat1s. These are approximations of the relative variances of hatZ_1^NldotshatZ_n^N = exp.(smcio.logZhats), i.e. the variances of hatZ_1^NhatZ_1ldotshatZ_n^NhatZ_n.When smcio.fullOutput = true, a vector of the approximations v^N_pn(f) or hatv^N_pn(f) for p in 1ldotsn can be obtained by callingSequentialMonteCarlo.vpns(smcio, f, hat::Bool, centred::Bool, n::Int64)One can choose n < smcio.n. If only the sum of one of these vectors is desired, i.e. v_n^N(f) or hatv^N_n(f), this can be obtained be callingSequentialMonteCarlo.v(smcio, f, hat::Bool, centred::Bool, n::Int64)"
},

{
    "location": "smcinterface.html#Adaptive-resampling-1",
    "page": "SMC interface",
    "title": "Adaptive resampling",
    "category": "section",
    "text": "The adaptive resampling mechanism is activated by choosing essThreshold <= 1.0. A vector of Bool values indicating whether or not resampling took place at each time can be accessed as smcio.resample, which has length smcio.n - 1 and is corresponds exactly to the random variables R_1 ldots R_n-1 defined in the SMC with adaptive resampling algorithm."
},

{
    "location": "smcinterface.html#Accessing-the-particle-system-1",
    "page": "SMC interface",
    "title": "Accessing the particle system",
    "category": "section",
    "text": "If smcio.fullOutput == true, one can access:name value\nsmcio.allZetas[p] zeta_p^1 ldots zeta_p^N\nsmcio.allWs[p] propto G_p(zeta_p^1) ldots G_p(zeta_p^N)\nsmcio.allEves[p] E_p^1 ldots E_p^N\nsmcio.allAs[p] A_p^1 ldots A_p^NNote that smcio.allAs has length smcio.n-1 while the others in the table above have length smcio.n.Even if smcio.fullOutput == false, one can access:name value\nsmcio.zetas zeta_n^1 ldots zeta_n^N\nsmcio.ws propto G_n(zeta_n^1)ldotsG_n(zeta_n^N)\nsmcio.eves E_n^1 ldots E_n^N\nsmcio.esses mathcalE_1^N ldots mathcalE_n^N"
},

{
    "location": "smcinterface.html#Conditional-SMC-1",
    "page": "SMC interface",
    "title": "Conditional SMC",
    "category": "section",
    "text": "The cSMC algorithm can be called as follows:csmc!(model::SMCModel, smcio::SMCIO, ref::Vector{Particle},\n  refout::Vector{Particle})where ref is the input reference path and refout the output path. It is permitted for ref and refout to be the same vector."
},

{
    "location": "guide.html#",
    "page": "Types and functions",
    "title": "Types and functions",
    "category": "page",
    "text": ""
},

{
    "location": "guide.html#Documentation-strings-1",
    "page": "Types and functions",
    "title": "Documentation strings",
    "category": "section",
    "text": ""
},

{
    "location": "guide.html#SequentialMonteCarlo.SMCModel",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.SMCModel",
    "category": "Type",
    "text": "SMCModel(M!::F1, lG::F2, maxn::Int64, particle::Type, pScratch::Type) where\n  {F1<:Function,F2<:Function}\n\nM! Mutation function\nlG Log potential function\nmaxn Maximum n for which the model is well-defined\nparticle Type of a particle\npScratch Type of particle scratch space\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.SMCIO",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.SMCIO",
    "category": "Type",
    "text": "SMCIO{Particle, ParticleScratch}\n\nStructs of this type should be constructed using the provided constructor. Important fields:\n\nN::Int64 Number of particles N\nn::Int64 Number of steps n\nnthreads::Int64 Number of threads\nfullOutput::Bool Whether particle system history should be recorded\nessThreshold::Float64 Relative ESS Threshold tau\nzetas::Vector{Particle} Time n particles zeta_n^1 ldots zeta_n^N\neves::Vector{Int64} Time n Eve indices E_n^1 ldots E_n^N\nws::Vector{Float64} Time n weights W_n^1 ldots W_n^N\nlogZhats::Vector{Float64} log(hatZ^N_1) ldots log(hatZ^N_n)\nVhat1s::Vector{Float64} hatV_1^N(1) ldots hatV_n^N(1)\nesses::Vector{Float64} Relative ESS values mathcalE_1^N ldots mathcalE_n^N\nresample::Vector{Bool} Resampling indicators R_1 ldots R_n-1\n\nPopulated only if fullOutput == true\n\nallZetas::Vector{Vector{Particle}} All the particles\nallWs::Vector{Vector{Float64}} All the weights\nallAs::Vector{Vector{Int64}} All the ancestor indices\nallEves::Vector{Vector{Int64}} All the Eve indices\n\n\n\n"
},

{
    "location": "guide.html#Types-1",
    "page": "Types and functions",
    "title": "Types",
    "category": "section",
    "text": "SMCModel{F1<:Function,F2<:Function}SMCIO{Particle, ParticleScratch}"
},

{
    "location": "guide.html#SequentialMonteCarlo.smc!-Tuple{SequentialMonteCarlo.SMCModel,SequentialMonteCarlo.SMCIO}",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.smc!",
    "category": "Method",
    "text": "smc!(model::SMCModel, smcio::SMCIO)\n\nRun the SMC algorithm for the given model and input/output arguments.\n\nIf smcio.nthreads == 1 the algorithm will run in serial.\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.csmc!-Union{Tuple{Particle}, Tuple{SequentialMonteCarlo.SMCModel,SequentialMonteCarlo.SMCIO{Particle,ParticleScratch} where ParticleScratch,Array{Particle,1},Array{Particle,1}}} where Particle",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.csmc!",
    "category": "Method",
    "text": "csmc!(model::SMCModel, smcio::SMCIO, ref::Vector{Particle}, refout::Vector{Particle})\n\nRun the conditional SMC algorithm for the given model, input/output arguments, reference path and output path.\n\nIt is permitted for ref and refout to be the same. If smcio.nthreads == 1 the algorithm will run in serial.\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.SMCIO-Union{Tuple{Int64,Int64,Int64,Bool,Float64}, Tuple{Int64,Int64,Int64,Bool}, Tuple{ParticleScratch}, Tuple{Particle}} where ParticleScratch where Particle",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.SMCIO",
    "category": "Method",
    "text": "SMCIO{Particle, ParticleScratch}(N::Int64, n::Int64, nthreads::Int64,\n  fullOutput::Bool, essThreshold::Float64 = 2.0) where\n  {Particle, ParticleScratch}\n\nConstructor for SMCIO structs.\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.eta-Union{Tuple{F}, Tuple{Particle}, Tuple{SequentialMonteCarlo.SMCIO{Particle,ParticleScratch} where ParticleScratch,F,Bool,Int64}} where F<:Function where Particle",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.eta",
    "category": "Method",
    "text": "eta(smcio::SMCIO{Particle}, f::F, hat::Bool, p::Int64) where {Particle, F<:Function}\n\nCompute:\n\n!hat: eta^N_p(f)\nhat:  hateta_p^N(f)\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.allEtas-Union{Tuple{F}, Tuple{SequentialMonteCarlo.SMCIO,F,Bool}} where F<:Function",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.allEtas",
    "category": "Method",
    "text": "allEtas(smcio::SMCIO, f::F, hat::Bool) where F<:Function\n\nCompute eta(smcio::SMCIO, f::F, hat::Bool, p) for p in {1, …, smcio.n}\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.slgamma-Union{Tuple{F}, Tuple{SequentialMonteCarlo.SMCIO,F,Bool,Int64}} where F<:Function",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.slgamma",
    "category": "Method",
    "text": "slgamma(smcio::SMCIO, f::F, hat::Bool, p::Int64) where {Particle, F<:Function}\n\nCompute:\n\n!hat: (eta^N_p(f) geq 0 log gamma^N_p(f))\nhat:  (hateta^N_p(f) geq 0 log hatgamma_p^N(f))\n\nThe result is returned as a Tuple{Bool, Float64}: the first component represents whether the returned value is non-negative, the second is the logarithm of the absolute value of the approximation.\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.allGammas-Union{Tuple{F}, Tuple{SequentialMonteCarlo.SMCIO,F,Bool}} where F<:Function",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.allGammas",
    "category": "Method",
    "text": "allGammas(smcio::SMCIO, f::F, hat::Bool) where F<:Function\n\nCompute slgamma(smcio::SMCIO, f::F, hat::Bool, p) for p in {1, …, smcio.n}\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.V-Union{Tuple{F}, Tuple{Particle}, Tuple{SequentialMonteCarlo.SMCIO{Particle,ParticleScratch} where ParticleScratch,F,Bool,Bool,Int64}} where F<:Function where Particle",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.V",
    "category": "Method",
    "text": "V(smcio::SMCIO{Particle}, f::F, hat::Bool, centred::Bool, p::Int64) where\n  {Particle, F<:Function}\n\nCompute:\n\n!hat & !centred: V^N_p(f)\n!hat & centred: V^N_p(f-eta_p^N(f))\nhat & !centred: hatV_p^N(f)\nhat & centred:  hatV_p^N(f-hateta_p^N(f))\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.vpns-Union{Tuple{F}, Tuple{SequentialMonteCarlo.SMCIO,F,Bool,Bool,Int64}} where F<:Function",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.vpns",
    "category": "Method",
    "text": "vpns(smcio::SMCIO, f::F, hat::Bool, centred::Bool, n::Int64) where F<:Function\n\nCompute a vector of the values of, for p = 1, …, n,\n\n!hat & !centred: v^N_pn(f)\n!hat & centred:  v^N_pn(f-eta_n^N(f))\nhat & !centred:  hatv_pn^N(f)\nhat & centred:   hatv_pn^N(f-hateta_n^N(f))\n\nNote: if essThreshold <= 1.0, and resampling did not occur at every time, the length of the output will be less than n.\n\n\n\n"
},

{
    "location": "guide.html#SequentialMonteCarlo.v-Union{Tuple{F}, Tuple{SequentialMonteCarlo.SMCIO,F,Bool,Bool,Int64}} where F<:Function",
    "page": "Types and functions",
    "title": "SequentialMonteCarlo.v",
    "category": "Method",
    "text": "v(smcio::SMCIO, f::F, hat::Bool, centred::Bool, n::Int64) where F<:Function\n\nCompute:\n\n!hat & !centred: v^N_n(f)\n!hat & centred: v^N_n(f-eta_n^N(f))\nhat & !centred: hatv_n^N(f)\nhat & centred:  hatv_n^N(f-hateta_n^N(f))\n\n\n\n"
},

{
    "location": "guide.html#Functions-1",
    "page": "Types and functions",
    "title": "Functions",
    "category": "section",
    "text": "smc!(model::SMCModel, smcio::SMCIO)csmc!(model::SMCModel, smcio::SMCIO{Particle}, ref::Vector{Particle},\n  refout::Vector{Particle}) where ParticleSMCIO{Particle, ParticleScratch}(N::Int64, n::Int64, nthreads::Int64,\n  fullOutput::Bool, essThreshold::Float64 = 2.0) where {Particle,\n  ParticleScratch}SequentialMonteCarlo.eta(smcio::SMCIO{Particle}, f::F, hat::Bool,\n  p::Int64) where {Particle, F<:Function}SequentialMonteCarlo.allEtas(smcio::SMCIO, f::F, hat::Bool) where F<:FunctionSequentialMonteCarlo.slgamma(smcio::SMCIO, f::F, hat::Bool, p::Int64) where\n  F<:FunctionSequentialMonteCarlo.allGammas(smcio::SMCIO, f::F, hat::Bool) where F<:FunctionSequentialMonteCarlo.V(smcio::SMCIO{Particle}, f::F, hat::Bool, centred::Bool,\n  p::Int64) where {Particle, F<:Function}SequentialMonteCarlo.vpns(smcio::SMCIO, f::F, hat::Bool, centred::Bool,\n  n::Int64) where F<:FunctionSequentialMonteCarlo.v(smcio::SMCIO, f::F, hat::Bool, centred::Bool,\n  n::Int64) where F<:Function"
},

{
    "location": "refs.html#",
    "page": "References",
    "title": "References",
    "category": "page",
    "text": ""
},

{
    "location": "refs.html#References-1",
    "page": "References",
    "title": "References",
    "category": "section",
    "text": "C. Andrieu, A. Doucet, and R. Holenstein. Particle Markov chain Monte Carlo methods. J. R. Stat. Soc. Ser. B Stat. Methodol., 72(3):269–342, 2010.\nC. Andrieu, A. Lee, and M. Vihola. Uniform ergodicity of the iterated conditional SMC and geometric ergodicity of particle Gibbs samplers. Bernoulli, 24(2):842–-872, 2018.\nH. P. Chan and T. L. Lai. A general theory of particle filters in hidden Markov models and some applications. Ann. Statist. 41(6):2877–2904, 2013.\nN. Chopin and S. S. Singh. On particle Gibbs sampling. Bernoulli, 21(3):1855–-1883, 2015.\nP. Del Moral. Feynman–Kac formulae: genealogical and interacting particle systems with applications. Springer-Verlag, 2004.\nP. Del Moral, A. Doucet, and A. Jasra. On adaptive resampling procedures for sequential Monte Carlo methods. Bernoulli, 18(1):252–278, 2012.\nL. Devroye. Non-uniform random variate generation. Springer-Verlag, 1986.\nA. Doucet and A. M. Johansen. A tutorial on particle filtering and smoothing: Fifteen years later. In D. Crisan and B. Rozovsky, eds., The Oxford Handbook of Nonlinear Filtering, pp. 656–-704. Oxford University Press, 2011.\nA. Doucet and A. Lee.  Sequential Monte Carlo methods. In M. Drton, S. Lauritzen, M. Maathuis, and M. Wainwright, eds., Handbook of Graphical Models. 2018. In press.\nN. J. Gordon, D. J. Salmond, and A. F. M. Smith. Novel approach to nonlinear/non-Gaussian Bayesian state estimation. Radar and Signal Processing, IEE Proceedings F, 140(2):107–-113, 1993.\nW. Hörmann. The generation of binomial random variates. J. Stat. Comput. Simul., 46(1–2):101–-110, 1993.\nG. Kitagawa. A Monte Carlo filtering and smoothing method for non-Gaussian nonlinear state space models. In Proceedings of the 2nd US-Japan Joint Seminar on Statistical Time Series Analysis, pp. 110–131, 1993.\nA. Kong, J. S. Liu, and W. H. Wong. Sequential imputations and Bayesian missing data problems. J. Am. Stat. Assoc., 89(425):278–288, 1994.\nA. Lee and N. Whiteley. Variance estimation in the particle filter. arXiv:1509.00394, 2015.\nA. Lee, C. Yau, M. B. Giles, A. Doucet, and C. C. Holmes. On the utility of graphics cards to perform massively parallel simulation of advanced Monte Carlo methods. J. Comput. Graph. Statist., 19(4):769–-789, 2010.\nF. Lindsten, R. Douc, and E. Moulines. Uniform ergodicity of the particle Gibbs sampler. Scand. J. Statist., 42(3):775-–797, 2015.\nJ. S. Liu and R. Chen. Blind deconvolution via sequential imputations. J. Am. Stat. Assoc., 90:567–-576, 1995.\nD. Lurie and H. O. Hartley. Machine-generation of order statistics for Monte Carlo computations. Am. Stat., 26(1):26–-27, 1972.\nL. M. Murray. Bayesian state-space modelling on high-performance hardware using LibBi. J. Stat. Softw. 67(10):1–36, 2015.\nL. Stewart and P. McCarty Jr. Use of Bayesian belief networks to fuse continuous and discrete information for target recognition, tracking, and situation assessment. In Aerospace Sensing, pp. 177-–185. International Society for Optics and Photonics, 1992.\nN. Whiteley, A. Lee, and K. Heine. On the role of interaction in sequential Monte Carlo algorithms. Bernoulli, 22(1):494–-529, 2016."
},

]}
