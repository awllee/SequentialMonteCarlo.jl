# SequentialMonteCarlo.jl Documentation

## Types

```@docs
SMCModel{F1<:Function,F2<:Function}
```

```@docs
SMCIO{Particle, ParticleScratch}
```

## Functions

```@docs
smc!(model::SMCModel, smcio::SMCIO)
```

```@docs
csmc!(model::SMCModel, smcio::SMCIO{Particle}, ref::Vector{Particle},
  refout::Vector{Particle}) where Particle
```

```@docs
SMCIO{Particle, ParticleScratch}(N::Int64, n::Int64, nthreads::Int64,
  fullOutput::Bool, essThreshold::Float64 = 2.0) where {Particle,
  ParticleScratch}
```

```@docs
SequentialMonteCarlo.eta(smcio::SMCIO{Particle}, f::F, hat::Bool,
  p::Int64) where {Particle, F<:Function}
```

```@docs
SequentialMonteCarlo.allEtas(smcio::SMCIO, f::F, hat::Bool) where F<:Function
```

```@docs
SequentialMonteCarlo.slgamma(smcio::SMCIO, f::F, hat::Bool, p::Int64) where
  F<:Function
```

```@docs
SequentialMonteCarlo.allGammas(smcio::SMCIO, f::F, hat::Bool) where F<:Function
```

```@docs
SequentialMonteCarlo.V(smcio::SMCIO{Particle}, f::F, hat::Bool, centred::Bool,
  p::Int64) where {Particle, F<:Function}
```

```@docs
SequentialMonteCarlo.vpns(smcio::SMCIO, f::F, hat::Bool, centred::Bool,
  n::Int64) where F<:Function
```

```@docs
SequentialMonteCarlo.v(smcio::SMCIO, f::F, hat::Bool, centred::Bool,
  n::Int64) where F<:Function
```
