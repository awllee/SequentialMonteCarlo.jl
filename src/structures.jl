import Compat.hasmethod

"""
    SMCModel(M!::F1, lG::F2, maxn::Int64, particle::Type, pScratch::Type) where
      {F1<:Function,F2<:Function}
- ```M!``` Mutation function
- ```lG``` Log potential function
- ```maxn``` Maximum n for which the model is well-defined
- ```particle``` Type of a particle
- ```pScratch``` Type of particle scratch space
"""
struct SMCModel{F1<:Function,F2<:Function}
  M!::F1
  lG::F2
  maxn::Int64
  particle::Type
  pScratch::Type
end

const SubVectorF64 = SubArray{Float64, 1, Array{Float64, 1},
  Tuple{UnitRange{Int64}}, true}
const SubVectorI64 = SubArray{Int64, 1, Array{Int64, 1},
  Tuple{UnitRange{Int64}}, true}
const SubVector{Particle} = SubArray{Particle,1,Array{Particle,1},Tuple{UnitRange{Int64}},true}

# the fields of this struct are used by the parallel SMC algorithm
# they are not intended for use by users of the package
mutable struct _SMCInternalParallel{Particle, ParticleScratch}
  Nperthread::Int64
  maxlws::Vector{Float64}
  partialSums::Vector{Float64}
  partialSumSqs::Vector{Float64}
  Ns::Vector{Int64}
  NsPartial::Vector{Int64}
  particleScratches::Vector{ParticleScratch}

  localAs::Vector{SubVectorI64}
  localEves::Vector{SubVectorI64}
  localOldEves::Vector{SubVectorI64}
  localZetas::Vector{SubVector{Particle}}
  localZetaAncs::Vector{SubVector{Particle}}
  localWs::Vector{SubVectorF64}
  localLogWs::Vector{SubVectorF64}
  localScratch1s::Vector{SubVectorF64}

  ## below are only populated if fullOutput = true
  localAllAs::Vector{Vector{SubVectorI64}}
  localAllZetas::Vector{Vector{SubVector{Particle}}}
  localAllEves::Vector{Vector{SubVectorI64}}
  localAllWs::Vector{Vector{SubVectorF64}}
end

function _SMCInternalParallel{Particle, ParticleScratch}(N::Int64, n::Int64,
  nthreads::Int64, fullOutput::Bool) where {Particle, ParticleScratch}
  @assert mod(N, nthreads) == 0 "N must be a multiple of nthreads"
  Nperthread = div(N, nthreads)

  if nthreads > 1
    maxlws::Vector{Float64} = Vector{Float64}(undef, nthreads)
    partialSums::Vector{Float64} = Vector{Float64}(undef, nthreads)
    partialSumSqs::Vector{Float64} = Vector{Float64}(undef, nthreads)
    Ns::Vector{Int64} = Vector{Int64}(undef, nthreads)
    NsPartial::Vector{Int64} = Vector{Int64}(undef, nthreads)

    particleScratches::Vector{ParticleScratch} =
      Vector{ParticleScratch}(undef, nthreads)
    # avoid false sharing
    Threads.@threads for i = 1:nthreads
      @inbounds particleScratches[i] = ParticleScratch()
    end

    localAs::Vector{SubVectorI64} =
      Vector{SubVectorI64}(undef, nthreads)
    localEves::Vector{SubVectorI64} =
      Vector{SubVectorI64}(undef, nthreads)
    localOldEves::Vector{SubVectorI64} =
      Vector{SubVectorI64}(undef, nthreads)
    localZetas::Vector{SubVector{Particle}} =
      Vector{SubVector{Particle}}(undef, nthreads)
    localZetaAncs::Vector{SubVector{Particle}} =
      Vector{SubVector{Particle}}(undef, nthreads)
    localWs::Vector{SubVectorF64} =
      Vector{SubVectorF64}(undef, nthreads)
    localLogWs::Vector{SubVectorF64} =
      Vector{SubVectorF64}(undef, nthreads)
    localScratch1s::Vector{SubVectorF64} =
      Vector{SubVectorF64}(undef, nthreads)
  else
    maxlws = Vector{Float64}(undef, 0)
    partialSums = Vector{Float64}(undef, 0)
    partialSumSqs = Vector{Float64}(undef, 0)
    Ns = Vector{Int64}(undef, 0)
    NsPartial = Vector{Int64}(undef, 0)
    particleScratches = Vector{ParticleScratch}(undef, 0)

    localAs = Vector{SubVectorI64}(undef, 0)
    localEves = Vector{SubVectorI64}(undef, 0)
    localOldEves = Vector{SubVectorI64}(undef, 0)
    localZetas = Vector{SubVector{Particle}}(undef, 0)
    localZetaAncs = Vector{SubVector{Particle}}(undef, 0)
    localWs = Vector{SubVectorF64}(undef, 0)
    localLogWs = Vector{SubVectorF64}(undef, 0)
    localScratch1s = Vector{SubVectorF64}(undef, 0)
  end

  if fullOutput && nthreads > 1
    localAllZetas::Vector{Vector{SubVector{Particle}}} =
      Vector{Vector{SubVector{Particle}}}(undef, n)
    localAllAs::Vector{Vector{SubVectorI64}} =
      Vector{Vector{SubVectorI64}}(undef, n-1)
    localAllEves::Vector{Vector{SubVectorI64}} =
      Vector{Vector{SubVectorI64}}(undef, n)
    localAllWs::Vector{Vector{SubVectorF64}} =
      Vector{Vector{SubVectorF64}}(undef, n)
    for p = 1:n
      localAllZetas[p] = Vector{SubVector{Particle}}(undef, nthreads)
      localAllEves[p] = Vector{SubVectorI64}(undef, nthreads)
      localAllWs[p] = Vector{SubVectorF64}(undef, nthreads)
      p < n && (localAllAs[p] = Vector{SubVectorI64}(undef, nthreads))
    end
  else
    localAllAs = Vector{Vector{SubVectorI64}}(undef, 0)
    localAllZetas = Vector{Vector{SubVector{Particle}}}(undef, 0)
    localAllEves = Vector{Vector{SubVectorI64}}(undef, 0)
    localAllWs = Vector{Vector{SubVectorF64}}(undef, 0)
  end

  return _SMCInternalParallel(Nperthread, maxlws, partialSums, partialSumSqs,
    Ns, NsPartial, particleScratches, localAs, localEves, localOldEves,
    localZetas, localZetaAncs, localWs, localLogWs, localScratch1s, localAllAs,
    localAllZetas, localAllEves, localAllWs)
end

# the fields of this struct are used by the SMC algorithm
# they are not intended for use by users of the package
mutable struct _SMCInternal{Particle, ParticleScratch}
  zetaAncs::Vector{Particle}
  oldEves::Vector{Int64}
  as::Vector{Int64}
  lws::Vector{Float64}
  scratch1::Vector{Float64}
  scratch2::Vector{Float64}
  vecOnes::Vector{Float64}
  nresamples::Int64

  sws::Float64
  mws::Float64
  maxlw::Float64
  particleScratch::ParticleScratch

  parallel::_SMCInternalParallel{Particle, ParticleScratch}
end

## constructor for _SMCInternal
function _SMCInternal{Particle, ParticleScratch}(N::Int64, n::Int64,
  nthreads::Int64, fullOutput::Bool) where {Particle, ParticleScratch}

  lws = Vector{Float64}(undef, N)
  as = Vector{Int64}(undef, N)
  fill!(as, 1)
  scratch1 = Vector{Float64}(undef, N)
  scratch2 = Vector{Float64}(undef, N)
  vecOnes = ones(N)
  nresamples::Int64 = 0
  zetaAncs = Vector{Particle}(undef, N)
  for i in 1:N
    zetaAncs[i] = Particle()
  end
  oldEves = Vector{Int64}(undef, N)

  ## assign user-defined particle scratch space
  @assert ParticleScratch == Nothing || !isbits(ParticleScratch)
  particleScratch = ParticleScratch()

  parallel::_SMCInternalParallel{Particle, ParticleScratch} =
    _SMCInternalParallel{Particle, ParticleScratch}(N, n, nthreads, fullOutput)

  return _SMCInternal(zetaAncs, oldEves, as, lws, scratch1, scratch2, vecOnes,
    nresamples, 0.0, 0.0, 0.0, particleScratch, parallel)
end

function _assignThreadViews(internal::_SMCInternal{Particle},
  zetas::Vector{Particle}, eves::Vector{Int64}, ws::Vector{Float64},
  allZetas::Vector{Vector{Particle}}, allWs::Vector{Vector{Float64}},
  allAs::Vector{Vector{Int64}}, allEves::Vector{Vector{Int64}}) where Particle

  zetaAncs::Vector{Particle} = internal.zetaAncs
  as::Vector{Int64} = internal.as
  oldEves::Vector{Int64} = internal.oldEves
  lws::Vector{Float64} = internal.lws
  scratch1::Vector{Float64} = internal.scratch1

  ip::_SMCInternalParallel = internal.parallel
  nthreads::Int64 = length(ip.localWs)
  Nperthread::Int64 = ip.Nperthread
  for i = 1:nthreads
    start = Nperthread * (i - 1) + 1
    finish = start + Nperthread - 1
    ip.localAs[i] = view(as, start:finish)
    ip.localEves[i] = view(eves, start:finish)
    ip.localOldEves[i] = view(oldEves, start:finish)
    ip.localZetas[i] = view(zetas, start:finish)
    ip.localZetaAncs[i] = view(zetaAncs, start:finish)
    ip.localWs[i] = view(ws, start:finish)
    ip.localLogWs[i] = view(lws, start:finish)
    ip.localScratch1s[i] = view(scratch1, start:finish)
  end

  fullOutput::Bool = length(ip.localAllWs) > 0
  if fullOutput && nthreads > 1
    n::Int64 = length(allWs)
    for p = 1:n
      for i = 1:nthreads
        start = Nperthread * (i - 1) + 1
        finish = start + Nperthread - 1
        ip.localAllZetas[p][i] = view(allZetas[p], start:finish)
        ip.localAllEves[p][i] = view(allEves[p], start:finish)
        ip.localAllWs[p][i] = view(allWs[p], start:finish)
        p < n && (ip.localAllAs[p][i] = view(allAs[p], start:finish))
      end
    end
  end
end

## SMC input / output struct. Also contains internal state for the SMC
## implementation
"""
    SMCIO{Particle, ParticleScratch}
Structs of this type should be constructed using the provided constructor.
Important fields:
- ```N::Int64``` Number of particles ``N``
- ```n::Int64``` Number of steps ``n``
- ```nthreads::Int64``` Number of threads
- ```fullOutput::Bool``` Whether particle system history should be recorded
- ```essThreshold::Float64``` Relative ESS Threshold ``\\tau``
- ```zetas::Vector{Particle}``` Time n particles ``\\zeta_n^1, \\ldots, \\zeta_n^N``
- ```eves::Vector{Int64}``` Time n Eve indices ``E_n^1, \\ldots, E_n^N``
- ```ws::Vector{Float64}``` Time n weights ``W_n^1, \\ldots, W_n^N``
- ```logZhats::Vector{Float64}``` ``\\log(\\hat{Z}^N_1), \\ldots, \\log(\\hat{Z}^N_n)``
- ```Vhat1s::Vector{Float64}``` ``\\hat{V}_1^N(1), \\ldots, \\hat{V}_n^N(1)``
- ```esses::Vector{Float64}``` Relative ESS values ``\\mathcal{E}_1^N, \\ldots, \\mathcal{E}_n^N``
- ```resample::Vector{Bool}``` Resampling indicators ``R_1, \\ldots, R_{n-1}``
Populated only if ```fullOutput == true```
- ```allZetas::Vector{Vector{Particle}}``` All the particles
- ```allWs::Vector{Vector{Float64}}``` All the weights
- ```allAs::Vector{Vector{Int64}}``` All the ancestor indices
- ```allEves::Vector{Vector{Int64}}``` All the Eve indices
"""
struct SMCIO{Particle, ParticleScratch}
  N::Int64
  n::Int64
  nthreads::Int64
  zetas::Vector{Particle}
  eves::Vector{Int64}
  ws::Vector{Float64}
  logZhats::Vector{Float64}
  Vhat1s::Vector{Float64}
  esses::Vector{Float64}
  resample::Vector{Bool}
  fullOutput::Bool
  essThreshold::Float64

  internal::_SMCInternal{Particle, ParticleScratch} # for internal use only

  # these are only populated if fullOutput = true
  allZetas::Vector{Vector{Particle}}
  allWs::Vector{Vector{Float64}}
  allAs::Vector{Vector{Int64}}
  allEves::Vector{Vector{Int64}}
end

"""
    SMCIO{Particle, ParticleScratch}(N::Int64, n::Int64, nthreads::Int64,
      fullOutput::Bool, essThreshold::Float64 = 2.0) where
      {Particle, ParticleScratch}
Constructor for ```SMCIO``` structs.
"""
function SMCIO{Particle, ParticleScratch}(N::Int64, n::Int64, nthreads::Int64,
  fullOutput::Bool, essThreshold::Float64 = 2.0) where {Particle,
  ParticleScratch}
  @assert hasmethod(Particle, ()) "Particle() must exist"
  @assert hasmethod(ParticleScratch, ()) "ParticleScratch() must exist"

  zetas = Vector{Particle}(undef, N)
  for i=1:N
    zetas[i] = Particle()
  end

  eves = Vector{Int64}(undef, N)
  ws = Vector{Float64}(undef, N)
  logZhats = Vector{Float64}(undef, n)
  Vhat1s = Vector{Float64}(undef, n)
  esses = Vector{Float64}(undef, n)
  resample = Vector{Bool}(undef, n-1)

  if fullOutput
    allZetas = Vector{Vector{Particle}}(undef, n)
    for i=1:n
      allZetas[i] = Vector{Particle}(undef, N)
      for j=1:N
        allZetas[i][j] = Particle()
      end
    end
    allWs = Vector{Vector{Float64}}(undef, n)
    for i=1:n
      allWs[i] = Vector{Float64}(undef, N)
    end
    allAs = Vector{Vector{Int64}}(undef, n-1)
    for i=1:n-1
      allAs[i] = Vector{Int64}(undef, N)
    end
    allEves = Vector{Vector{Int64}}(undef, n)
    for i=1:n
      allEves[i] = Vector{Int64}(undef, N)
    end
  else
    allZetas = Vector{Vector{Particle}}(undef, 0)
    allWs = Vector{Vector{Float64}}(undef, 0)
    allAs = Vector{Vector{Int64}}(undef, 0)
    allEves = Vector{Vector{Int64}}(undef, 0)
  end

  internal::_SMCInternal = _SMCInternal{Particle, ParticleScratch}(N, n,
    nthreads, fullOutput)

  _assignThreadViews(internal, zetas, eves, ws, allZetas, allWs, allAs, allEves)

  return SMCIO(N, n, nthreads, zetas, eves, ws, logZhats, Vhat1s, esses,
    resample, fullOutput, essThreshold, internal, allZetas, allWs, allAs,
    allEves)
end
