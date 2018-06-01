import NonUniformRandomVariateGeneration.sampleMultinomial!

@inline function _iotaParallel!(array::Vector{Int64}, N::Int64, nthreads::Int64,
  Nperthread::Int64)
  Threads.@threads for i in 1:nthreads
    start::Int64 = Nperthread * (i - 1) + 1
    finish::Int64 = start + Nperthread - 1
    for i = start:finish
      @inbounds array[i] = i
    end
  end
end

@inline function _cpeHelper(smcio::SMCIO)
  Threads.@threads for i = 1:smcio.nthreads
    ip::_SMCInternalParallel = smcio.internal.parallel
    oldEves::Vector{Int64} = smcio.internal.oldEves
    @inbounds localAs::SubVectorI64 = ip.localAs[i]
    @inbounds localEves::SubVectorI64 = ip.localEves[i]
    _setEves!(localEves, oldEves, localAs)
  end
end

@inline function _copyParticlesEves!(smcio::SMCIO{Particle}) where Particle
  Threads.@threads for i = 1:smcio.nthreads
    ip::_SMCInternalParallel = smcio.internal.parallel
    zetas::Vector{Particle} = smcio.zetas
    @inbounds localZetaAncs::SubVector{Particle} = ip.localZetaAncs[i]
    @inbounds localAs::SubVectorI64 = ip.localAs[i]
    @inbounds localOldEves::SubVectorI64 = ip.localOldEves[i]
    @inbounds localEves::SubVectorI64 = ip.localEves[i]
    _copyParticles!(localZetaAncs, zetas, localAs)
    localOldEves .= localEves
  end
  _cpeHelper(smcio)
end

@inline function _MutateWeightHelper1(p::Int64, smcio::SMCIO)
  ip::_SMCInternalParallel= smcio.internal.parallel
  maxlw::Float64 = maximum(ip.maxlws)
  smcio.internal.maxlw = maxlw

  for i = 1:smcio.nthreads
    @inbounds tmp::Float64 = exp(ip.maxlws[i] - maxlw)
    @inbounds ip.partialSums[i] *= tmp
    @inbounds ip.partialSumSqs[i] *= tmp * tmp
  end

  smcio.internal.sws = sum(ip.partialSums)
  smcio.internal.mws = smcio.internal.sws / smcio.N
  @inbounds smcio.esses[p] = smcio.internal.mws * smcio.internal.mws *
    smcio.N / sum(ip.partialSumSqs)
end

@inline function _MutateWeightHelper2(smcio::SMCIO)
  Threads.@threads for i = 1:smcio.nthreads
    ip::_SMCInternalParallel = smcio.internal.parallel
    @inbounds G::Float64 = exp(ip.maxlws[i] - smcio.internal.maxlw) /
      smcio.internal.mws
    @inbounds ip.localWs[i] .*= G
  end
end

@inline function _MutateWeight!(p::Int64, model::SMCModel,
  smcio::SMCIO{Particle, ParticleScratch},
  ref::Union{Nothing, Particle} = nothing) where {Particle, ParticleScratch}
  Threads.@threads for i = 1:smcio.nthreads
    ip::_SMCInternalParallel{Particle, ParticleScratch} = smcio.internal.parallel
    @inbounds localZetas::SubVector{Particle} = ip.localZetas[i]
    @inbounds localZetaAncs::SubVector{Particle} = ip.localZetaAncs[i]
    @inbounds localLogWs::SubVectorF64 = ip.localLogWs[i]
    @inbounds localWs::SubVectorF64 = ip.localWs[i]
    @inbounds pScratch::ParticleScratch = ip.particleScratches[i]
    engine::SMCRNG = getSMCRNG()

    if i == 1 && ref != nothing
      _mutateParticles!(localZetas, engine, p, model.M!, localZetaAncs,
        pScratch, ref)
    else
      _mutateParticles!(localZetas, engine, p, model.M!, localZetaAncs,
        pScratch)
    end
    @inbounds noresample::Bool = p == 1 || smcio.resample[p-1]
    _logWeightParticles!(localLogWs, p, model.lG, localZetas, pScratch,
      noresample, localWs)

    localmaxlw::Float64 = maximum(localLogWs)
    @inbounds ip.maxlws[i] = localmaxlw
    localWs .= exp.(localLogWs .- localmaxlw)
    @inbounds ip.partialSums[i] = sum(localWs)
    @inbounds ip.partialSumSqs[i] = mapreduce(x -> x * x, +, localWs)
  end

  _MutateWeightHelper1(p, smcio)
  _MutateWeightHelper2(smcio)
end

@inline function _resampleHelper(smcio::SMCIO, conditional::Bool)
  Threads.@threads for i = 1:smcio.nthreads
    ip::_SMCInternalParallel = smcio.internal.parallel
    @inbounds if ip.Ns[i] > 0
      @inbounds localScratch1 = ip.localScratch1s[i]
      @inbounds localWs = ip.localWs[i]
      @inbounds localN = ip.Ns[i]
      engine::SMCRNG = getSMCRNG()
      start::Int64 = ip.Nperthread * (i - 1)
      @inbounds offset::Int64 = i == 1 ? 0 : ip.NsPartial[i-1]
      if i == 1 && conditional
        @inbounds smcio.internal.as[1] = 1
        if localN > 1
          sampleSerial!(smcio.internal.as, engine, localWs, localN - 1,
            smcio.internal.scratch2, localScratch1, offset + 1, start)
        end
      else
        sampleSerial!(smcio.internal.as, engine, localWs, localN,
          smcio.internal.scratch2, localScratch1, offset, start)
      end
    end
  end
end

@inline function _resampleParallel!(smcio::SMCIO, p::Int64, conditional::Bool)
  if p != smcio.n-1 && smcio.esses[p] > smcio.essThreshold
    smcio.resample[p] = false
    _iotaParallel!(smcio.internal.as, smcio.N, smcio.nthreads,
      smcio.internal.parallel.Nperthread)
    return
  end
  @inbounds smcio.resample[p] = true
  smcio.internal.nresamples += 1

  ip::_SMCInternalParallel = smcio.internal.parallel
  Ns::Vector{Int64} = ip.Ns
  NsPartial::Vector{Int64} = ip.NsPartial
  partialSums::Vector{Float64} = ip.partialSums
  partialSums ./= smcio.internal.sws
  if conditional
    sampleMultinomial!(smcio.N-1, partialSums, Ns, getSMCRNG())
    @inbounds Ns[1] += 1
  else
    sampleMultinomial!(smcio.N, partialSums, Ns, getSMCRNG())
  end

  cumsum!(NsPartial,Ns)
  _resampleHelper(smcio, conditional)
end

@inline function _Vhat1ParallelHelper(smcio::SMCIO)
  N::Int64 = smcio.N
  vs::Vector{Float64} = smcio.internal.scratch1

  v::Float64 = 0.0
  for t = 1:smcio.nthreads
    @inbounds v += vs[t]
  end
  q::Int64 = 1 + smcio.internal.nresamples
  return 1 - (1-v)*(N/(N-1))^q
end

@inline function _Vhat1Parallel(smcio::SMCIO)
  # N::Int64 = smcio.N
  # vs::Vector{Float64} = smcio.internal.scratch1
  Threads.@threads for t = 1:smcio.nthreads
    N::Int64 = smcio.N
    vs::Vector{Float64} = smcio.internal.scratch1
    Nperthread::Int64 = smcio.internal.parallel.Nperthread
    eves::Vector{Int64} = smcio.eves
    ws::Vector{Float64} = smcio.ws
    i::Int64 = searchsortedfirst(eves, Nperthread * (t - 1) + 1)
    finish::Int64 = searchsortedlast(eves, Nperthread * t)
    vLocal::Float64 = 0.0
    while i <= finish
      @inbounds e = eves[i]
      tmp::Float64 = 0.0
      @inbounds while (i <= N && eves[i] == e)
        @inbounds tmp += ws[i]
        i += 1
      end
      tmp /= N
      vLocal += tmp*tmp
    end
    vs[t] = vLocal
  end
  return _Vhat1ParallelHelper(smcio)
end

@inline function _intermediateOutputParallel!(smcio::SMCIO{Particle},
  p::Int64) where Particle
  Threads.@threads for i = 1:smcio.nthreads
    ip::_SMCInternalParallel{Particle} = smcio.internal.parallel
    @inbounds _copyParticles!(ip.localAllZetas[p][i], ip.localZetas[i])
    @inbounds ip.localAllWs[p][i] .= ip.localWs[i]
    @inbounds ip.localAllEves[p][i] .= ip.localEves[i]
    @inbounds p < smcio.n && (ip.localAllAs[p][i] .= ip.localAs[i])
  end
end

# Parallel implementation of SMC
@inline function _smcParallel!(model::SMCModel, smcio::SMCIO{Particle},
  ref::Union{Nothing, Vector{Particle}} = nothing,
  refout::Union{Nothing, Vector{Particle}} = nothing) where Particle
  smcio.internal.nresamples = 0
  lZ::Float64 = 0.0
  _iotaParallel!(smcio.eves, smcio.N, smcio.nthreads,
    smcio.internal.parallel.Nperthread)

  for p = 1:smcio.n
    p > 1 && _copyParticlesEves!(smcio)
    if ref == nothing
      _MutateWeight!(p, model, smcio)
    else
      _MutateWeight!(p, model, smcio, ref[p])
    end

    lZ += smcio.internal.maxlw + log(smcio.internal.mws)
    @inbounds smcio.logZhats[p] = lZ

    @inbounds smcio.Vhat1s[p] = _Vhat1Parallel(smcio)

    p < smcio.n && _resampleParallel!(smcio, p, ref != nothing)
    smcio.fullOutput && _intermediateOutputParallel!(smcio, p)
  end

  refout != nothing && SequentialMonteCarlo.pickParticle!(refout, smcio)
  return
end
