@inline function _iotaParallel!(array::Vector{Int64})
  Threads.@threads for i in 1:length(array)
    @inbounds array[i] = i
  end
end

@inline function _cpeHelper(smcio::SMCIO)
  Threads.@threads for i = 1:smcio.nthreads
    @inbounds localAs::SubVectorI64 = smcio.internal.localAs[i]
    @inbounds localEves::SubVectorI64 = smcio.internal.localEves[i]
    _setEves!(localEves, smcio.internal.oldEves, localAs)
  end
end

@inline function _copyParticlesEves!(model::SMCModel, smcio::SMCIO)
  Threads.@threads for i = 1:smcio.nthreads
    @inbounds localZetaAncs = smcio.internal.localZetaAncs[i]
    @inbounds localAs = smcio.internal.localAs[i]
    @inbounds localOldEves = smcio.internal.localOldEves[i]
    @inbounds localEves = smcio.internal.localEves[i]
    _copyParticles!(localZetaAncs, smcio.zetas, localAs)
    @inbounds localOldEves .= localEves
  end
  _cpeHelper(smcio)
end

@inline function _MutateWeightHelper(smcio::SMCIO)
  Threads.@threads for i = 1:smcio.nthreads
    @inbounds localWs = smcio.internal.localWs[i]
    @inbounds G::Float64 = exp(smcio.internal.maxlws[i] -
      smcio.internal.maxlw) / smcio.internal.mws
    localWs .*= G
  end
end

@inline function _MutateWeight!(p::Int64, model::SMCModel,
  smcio::SMCIO{Particle}, ref::Union{Void, Particle} = nothing) where Particle
  nthreads::Int64 = smcio.nthreads
  pScratches::Vector{model.pScratch} = smcio.internal.particleScratches

  Threads.@threads for i = 1:nthreads
    @inbounds localZetas = smcio.internal.localZetas[i]
    @inbounds localZetaAncs = smcio.internal.localZetaAncs[i]
    @inbounds localLogWs = smcio.internal.localLogWs[i]
    @inbounds localWs = smcio.internal.localWs[i]
    engine::SMCRNG = getSMCRNG()

    if i == 1 && ref != nothing
      @inbounds _mutateParticles!(localZetas, engine, p, model.M!,
        localZetaAncs, pScratches[i], ref)
    else
      @inbounds _mutateParticles!(localZetas, engine, p, model.M!,
        localZetaAncs, pScratches[i])
    end
    @inbounds _logWeightParticles!(localLogWs, p, model.lG, localZetas,
      pScratches[i], p == 1 || smcio.resample[p-1], localWs)

    localmaxlw::Float64 = maximum(localLogWs)
    @inbounds smcio.internal.maxlws[i] = localmaxlw
    @inbounds localWs .= exp.(localLogWs .- localmaxlw)
    @inbounds smcio.internal.partialSums[i] = sum(localWs)
    @inbounds smcio.internal.partialSumSqs[i] = mapreduce(x -> x*x, +, localWs)
  end

  maxlw::Float64 = maximum(smcio.internal.maxlws)
  smcio.internal.maxlw = maxlw

  for i = 1:nthreads
    @inbounds tmp::Float64 = exp(smcio.internal.maxlws[i] - maxlw)
    @inbounds smcio.internal.partialSums[i] *= tmp
    @inbounds smcio.internal.partialSumSqs[i] *= tmp * tmp
  end

  smcio.internal.sws = sum(smcio.internal.partialSums)
  smcio.internal.mws = smcio.internal.sws / smcio.N
  @inbounds smcio.esses[p] = smcio.internal.mws * smcio.internal.mws *
    smcio.N / sum(smcio.internal.partialSumSqs)

  _MutateWeightHelper(smcio)
end

@inline function _resampleHelper(smcio::SMCIO, conditional::Bool)
  Ns::Vector{Int64} = smcio.internal.Ns
  NsPartial::Vector{Int64} = smcio.internal.NsPartial
  Nperthread::Int64 = smcio.internal.Nperthread
  Threads.@threads for i = 1:smcio.nthreads
    @inbounds if Ns[i] > 0
      @inbounds localScratch1 = smcio.internal.localScratch1s[i]
      @inbounds localWs = smcio.internal.localWs[i]
      engine::SMCRNG = getSMCRNG()
      start::Int64 = Nperthread * (i - 1)
      @inbounds offset::Int64 = i == 1 ? 0 : NsPartial[i-1]
      if i == 1 && conditional
        @inbounds smcio.internal.as[1] = 1
        @inbounds if Ns[i] > 1
          @inbounds sampleSerial!(smcio.internal.as, engine, localWs, Ns[i]-1,
            smcio.internal.scratch2, localScratch1, offset + 1, start)
        end
      else
        @inbounds sampleSerial!(smcio.internal.as, engine, localWs, Ns[i],
          smcio.internal.scratch2, localScratch1, offset, start)
      end
    end
  end
end

@inline function _resampleParallel!(smcio::SMCIO, p::Int64, conditional::Bool)
  if p != smcio.n-1 && smcio.esses[p] > smcio.essThreshold
    smcio.resample[p] = false
    _iotaParallel!(smcio.internal.as)
    return
  end
  @inbounds smcio.resample[p] = true
  smcio.internal.nresamples += 1

  Ns = smcio.internal.Ns
  NsPartial = smcio.internal.NsPartial
  partialSums = smcio.internal.partialSums
  partialSums ./= smcio.internal.sws
  if conditional
    @inbounds sampleMultinomial!(smcio.N-1, partialSums, Ns, getSMCRNG())
    @inbounds Ns[1] += 1
  else
    @inbounds sampleMultinomial!(smcio.N, partialSums, Ns, getSMCRNG())
  end

  cumsum!(NsPartial,Ns)
  _resampleHelper(smcio, conditional)
end

@inline function _Vhat1Parallel(smcio::SMCIO)
  N::Int64 = smcio.N
  eves::Vector{Int64} = smcio.eves
  ws::Vector{Float64} = smcio.ws
  q::Int64 = 1 + smcio.internal.nresamples

  nthreads::Int64 = smcio.nthreads
  Nperthread::Int64 = smcio.internal.Nperthread
  vs = smcio.internal.scratch1
  Threads.@threads for t = 1:nthreads
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
  v::Float64 = 0.0
  for t = 1:nthreads
    @inbounds v += vs[t]
  end
  return 1 - (1-v)*(N/(N-1))^q
end

@inline function _intermediateOutputParallel!(smcio::SMCIO{Particle},
  p::Int64) where Particle
  Threads.@threads for i = 1:smcio.nthreads
    _copyParticles!(smcio.internal.localAllZetas[p][i], smcio.internal.localZetas[i])
    @inbounds smcio.internal.localAllWs[p][i] .= smcio.internal.localWs[i]
    @inbounds smcio.internal.localAllEves[p][i] .= smcio.internal.localEves[i]
    if p < smcio.n
      @inbounds smcio.internal.localAllAs[p][i] .= smcio.internal.localAs[i]
    end
  end
end

# Parallel implementation of SMC
@inline function _smcParallel!(model::SMCModel, smcio::SMCIO{Particle},
  ref::Union{Void, Vector{Particle}} = nothing,
  refout::Union{Void, Vector{Particle}} = nothing) where Particle
  smcio.internal.nresamples = 0
  lZ::Float64 = 0.0
  _iotaParallel!(smcio.eves)

  for p = 1:smcio.n
    p > 1 && _copyParticlesEves!(model, smcio)
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
