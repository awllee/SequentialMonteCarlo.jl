# [Interface](@id interface)

## Specifying an SMC model

We recall from the [Introduction](@ref) that the SMC algorithm is defined in terms of $M_1, \ldots, M_n$ and $G_1, \ldots, G_n$.

A ```model``` of type ```SMCModel``` can be created by calling
```
model = SMCModel(M!, lG, maxn::Int64, Particle, ParticleScratch)
```
where ```Particle``` and ```ParticleScratch``` are user-defined types, ```M!``` is a void function with arguments
```
M!(newParticle::Particle, rng::SMCRNG, p::Int64, particle::Particle,
  scratch::ParticleScratch)
```
and ```lG``` is a function returning a ```Float64``` and has arguments
```
lG(p::Int64, particle::Particle, scratch::ParticleScratch)
```

There is a correspondence between the function ```M!``` and $M_1, \ldots, M_n$: calling ```M!(x', rng, p, x, scratch)```, should make $x'$ a realization of a sample from $M_p(x, \cdot)$ with the convention that $M_1(x,\cdot) = M_1(\cdot)$ for any $x$. Similarly, ```lG``` and $G_1, \ldots, G_n$ correspond in that ```lG(p,x)``` $ = \log G_p(x)$. Logarithms are used to avoid numerical issues. ```maxn``` is the maximum value of ```n``` for which ```M!``` and ```lG``` are well-defined; users may choose to run the SMC algorithm for any integer value of ```n``` up to and including ```maxn```.

The types ```Particle``` and ```ParticleScratch``` must have constructors that take no arguments. One may choose ```ParticleScratch = Void```, in which case ```nothing``` will be passed to ```M!``` and ```lG```. Using scratch space is optional but can significantly improve performance in certain scenarios; it provides a mechanism for users to avoid dynamic memory allocations in ```M!``` and/or ```lG```. This scratch space will be used by every particle associated with a given thread. A thread-specific pseudo-random number generator (RNG) ```rng``` will be passed to the ```M!``` function by the algorithm, and should be used in lieu of Julia's global RNG.

## Running the SMC algorithm

The SMC algorithm is run by calling
```
smc!(model::SMCModel, smcio::SMCIO)
```
The second argument, ```smcio```, is a struct containing inputs and outputs for the ```smc!``` algorithm. It is straightforward to construct and involves specifying the number of particles ```N```, the number of iterations ```n```, the number of threads ```nthreads```, whether the entire history of the particle system should be recorded ```fullOutput```, and an effective sample size threshold ```essThreshold``` that is explained on the [adaptive resampling page](@ref adaptiveresampling), and can be safely ignored on a first reading. ```smcio``` can be created by calling
```
smcio = SMCIO{model.particle, model.pScratch}(N::Int64, n::Int64,
  nthreads::Int64, fullOutput::Bool, essThreshold::Float64 = 2.0)
```
Currently ```N``` must be an integer multiple of ```nthreads```; this may change in the future.

## Extracting the main outputs

A vector of approximations $(\log\hat{Z}_1^N, \ldots, \log\hat{Z}_n^N)$ is stored in ```smcio.logZhats```, which can be used in conjunction with ```SequentialMonteCarlo.eta``` to produce the approximations $\gamma_p^N(f)$ or $\hat{\gamma}_p^N(f)$.

One can extract the approximation $\eta_p^N(f)$ or $\hat{\eta}_p^N(f)$ by calling
```
SequentialMonteCarlo.eta(smcio, f, hat::Bool, p::Int64)
```
with ```hat``` determining which approximation is returned. Calling this function with ```p < smcio.n``` requires ```smcio.fullOutput = true```.

One can extract the approximations $(\eta^N_p(f) \geq 0, \log |\gamma^N_p(f)|)$ or $(\hat{\eta}^N_p(f) \geq 0, \log |\hat{\gamma}^N_p(f)|)$ by calling
```
slgamma(smcio, f, hat::Bool, p::Int64)
```
with ```hat``` determining which approximation is returned. Calling this function with ```p < smcio.n``` requires ```smcio.fullOutput = true```.

## Variance estimators

The following functions relate to the approximations detailed in the [variance estimators page](@ref vepage).

The function
```
SequentialMonteCarlo.V(smcio, f, hat::Bool, centred::Bool, p::Int64)
```
returns:

| ```hat``` | ```centred``` | approximation | approximation of |
| --------- | ------------- | ------------- | ---------------- |
| false | false | $V_p^N(f)$ | ``{\rm var} \left \{ \gamma_p^N(f)/\gamma_p(1) \right \}`` |
| false | true  | $V_p^N(f-\eta_p^N(f))$ | ``\mathbb{E}\left[\left\{ \eta_p^N(f)-\eta_p(f)\right\} ^2\right]`` |
| true  | false | $\hat{V}_p^N(f)$ | ``{\rm var}\left\{ \hat{\gamma}_p^N(f)/\hat{\gamma}_p(1)\right\} `` |
| true  | true  | $\hat{V}_p^N(f-\hat{\eta}_p^N(f))$ | ``\mathbb{E}\left[\left\{ \hat{\eta}_p^N(f)-\hat{\eta}_p(f)\right\} ^2\right]`` |

Calling ```SequentialMonteCarlo.V``` with ```p < smcio.n``` requires ```smcio.fullOutput = true```.

A vector of the quantities $\hat{V}_1^N(1),\ldots,\hat{V}_n^N(1)$ is stored in ```smcio.Vhat1s```. These are approximations of the *relative* variances of $\hat{Z}_1^N,\ldots,\hat{Z}_n^N$ = ```exp.(smcio.logZhats)```, i.e. the variances of $\hat{Z}_1^N/\hat{Z}_1,\ldots,\hat{Z}_n^N/\hat{Z}_n$.

When ```smcio.fullOutput = true```, a vector of the approximations $v^N_{p,n}(f)$ or $\hat{v}^N_{p,n}(f)$ for $p \in \{1,\ldots,n\}$ can be obtained by calling
```
SequentialMonteCarlo.vpns(smcio, f, hat::Bool, centred::Bool, n::Int64)
```
One can choose ```n < smcio.n```. If only the sum of one of these vectors is desired, i.e. $v_n^N(f)$ or $\hat{v}^N_n(f)$, this can be obtained be calling
```
SequentialMonteCarlo.v(smcio, f, hat::Bool, centred::Bool, n::Int64)
```

## Adaptive resampling

The adaptive resampling mechanism is activated by choosing ```essThreshold <= 1.0```. A vector of ```Bool``` values indicating whether or not resampling took place at each time can be accessed as ```smcio.resample```, which has length ```smcio.n - 1``` and is corresponds exactly to the random variables $R_1, \ldots, R_{n-1}$ defined in the [SMC with adaptive resampling algorithm](@ref adaptiveresampling).

## Accessing the particle system

If ```smcio.fullOutput == true```, one can access:

| name                    | value |
| ----------------------- | ---------------- |
| ```smcio.allZetas[p]``` | ``\zeta_p^1, \ldots, \zeta_p^N`` |
| ```smcio.allWs[p]```    | ``\propto G_p(\zeta_p^1), \ldots, G_p(\zeta_p^N)`` |
| ```smcio.allEves[p]```  | ``E_p^1, \ldots, E_p^N`` |
| ```smcio.allAs[p]```    | ``A_p^1, \ldots, A_p^N`` |

Note that ```smcio.allAs``` has length ```smcio.n-1``` while the others in the table above have length ```smcio.n```.

Even if ```smcio.fullOutput == false```, one can access:

| name              | value |
| ----------------- | ---------------- |
| ```smcio.zetas``` | ``\zeta_n^1, \ldots, \zeta_n^N`` |
| ```smcio.ws```    | ``\propto G_n(\zeta_n^1),\ldots,G_n(\zeta_n^N)`` |
| ```smcio.eves```  | ``E_n^1, \ldots, E_n^N`` |
| ```smcio.esses```  | ``\mathcal{E}_{1}^N, \ldots, \mathcal{E}_{n}^N`` |

## Conditional SMC

The cSMC algorithm can be called as follows:
```
csmc!(model::SMCModel, smcio::SMCIO, ref::Vector{Particle},
  refout::Vector{Particle})
```
where ```ref``` is the input reference path and ```refout``` the output path. It is permitted for ```ref``` and ```refout``` to be the same vector.
