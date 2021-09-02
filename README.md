# SequentialMonteCarlo.jl

<!-- badges: start -->
[![CI](https://github.com/awllee/SequentialMonteCarlo.jl/workflows/CI/badge.svg)](https://github.com/awllee/SequentialMonteCarlo.jl/actions)
[![codecov](https://codecov.io/gh/awllee/SequentialMonteCarlo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/awllee/SequentialMonteCarlo.jl)
<!-- badges: end -->

This package provides a light interface to a serial and multi-threaded implementation of the Sequential Monte Carlo (SMC) algorithm. SMC is a random algorithm for approximate numerical integration and/or sampling.

The [documentation](https://awllee.github.io/SequentialMonteCarlo.jl/latest) and some [examples](https://github.com/awllee/SMCExamples.jl) may be helpful for getting started.

## Getting the package

```julia
Pkg.add("SequentialMonteCarlo")
```

## Quick start:

```julia
## Load the package

using SequentialMonteCarlo

## Define a particle consisting of one Float64

mutable struct Float64Particle
  x::Float64
  Float64Particle() = new()
end

## The initial distribution is a standard normal, and the Markov transitions
## define a non-stationary Markov chain, since 1.5 > 1.

function M!(newParticle::Float64Particle, rng, p::Int64,
  particle::Float64Particle, ::Nothing)
  if p == 1
    newParticle.x = randn(rng)
  else
    newParticle.x = 1.5 * particle.x + 0.5 * randn(rng)
  end
end

## The log potential function is x -> -x^2 so the potential function is
## x -> exp(-x^2).

function lG(p::Int64, particle::Float64Particle, ::Nothing)
  return - particle.x * particle.x
end

## This is a pedagogical example: one can deduce that eta_p is N(0, 1) for all
## p and \hat{eta}_p is N(0, 1/3) for all p. Essentially, the potential
## functions stop the Markov chain from exploding by favouring values
## closer to 0. In addition, \hat{Z}_p = (sqrt(3)/3)^p.

## Specify the model using M! and lG, stating that the maximum
## length of the model is 100, and specifying the types of the particle and
## particle scratch space. The latter is Nothing in this case as no scratch space
## is required.

model = SMCModel(M!, lG, 100, Float64Particle, Nothing)

## Create the SMC input/output struct, specifying the number of particles N as
## 2^20 = 1048576, that the algorithm should be run for 10 steps, that 1 thread
## should be used (i.e. it should run in serial) and the whole particle system
## should be recorded.

smcio = SMCIO{model.particle, model.pScratch}(2^20, 10, 1, true)

## Run the algorithm twice, timing it both times. The first time will include
## compilation overhead. The second time there will be no allocations (apart
## from the small number associated with using the @time macro).

@time smc!(model, smcio)
@time smc!(model, smcio)

## Check that the approximations in smcio.logZhats are close to the true values.

println(Vector(log(sqrt(3)/3) * (1:10)))
println(smcio.logZhats)

## Check that the first and second moments of the eta_p's (resp. \hat{eta}_p's)
## are close to 0 and 1 (resp. 0 and 1/3).

println(SequentialMonteCarlo.allEtas(smcio, p->p.x, false))
println(SequentialMonteCarlo.allEtas(smcio, p->p.x*p.x, false))

println(SequentialMonteCarlo.allEtas(smcio, p->p.x, true))
println(SequentialMonteCarlo.allEtas(smcio, p->p.x*p.x, true))

## Now try with Threads.nthreads() threads instead of 1 and time the algorithm.
## There are compilation overheads the first time for the parallel parts of the
## SMC implementation. There are still a number of small allocations the second
## time; there is an allocation for each parallel region in the algorithm, and
## there are a few such regions at each step of the algorithm. This is due to
## Julia's multi-threading interface.

nthreads = Threads.nthreads()
smcio = SMCIO{model.particle, model.pScratch}(2^20, 10, nthreads, true)

@time smc!(model, smcio)
@time smc!(model, smcio)
```
