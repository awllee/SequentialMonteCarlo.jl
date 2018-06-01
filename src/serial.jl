import Compat.Nothing

@inline function _iota!(array::Vector{Int64})
  for i in eachindex(array)
    @inbounds array[i] = i
  end
end

@inline function _intermediateOutputSerial!(smcio::SMCIO, p::Int64)
  if length(smcio.allZetas) > 0
    @inbounds vec = smcio.allZetas[p]
    _copyParticles!(vec, smcio.zetas)
    @inbounds smcio.allWs[p] .= smcio.ws
    @inbounds smcio.allEves[p] .= smcio.eves
    @inbounds p != smcio.n && (smcio.allAs[p] .= smcio.internal.as)
  end
end

function _Vhat1Serial(smcio::SMCIO)
  N::Int64 = smcio.N
  eves::Vector{Int64} = smcio.eves
  ws::Vector{Float64} = smcio.ws
  q::Int64 = 1 + smcio.internal.nresamples
  i::Int64 = 1
  v::Float64 = 0.0
  while i <= N
    @inbounds e = eves[i]
    tmp::Float64 = 0.0
    @inbounds while i <= N && eves[i] == e
      @inbounds tmp += ws[i]
      i += 1
    end
    tmp /= N
    v += tmp*tmp
  end
  return 1 - (1-v)*(N/(N-1))^q
end

@inline function _resampleSerial!(smcio::SMCIO, p::Int64, conditional::Bool)
  if p != smcio.n-1 && smcio.esses[p] > smcio.essThreshold
    smcio.resample[p] = false
    _iota!(smcio.internal.as)
    return
  end
  smcio.resample[p] = true
  smcio.internal.nresamples += 1
  if conditional
    smcio.internal.as[1] = 1
    sampleSerial!(smcio.internal.as, getSMCRNG(), smcio.ws, smcio.N-1,
      smcio.internal.scratch1, smcio.internal.scratch2, 1)
  else
    sampleSerial!(smcio.internal.as, getSMCRNG(), smcio.ws, smcio.N,
      smcio.internal.scratch1, smcio.internal.scratch2)
  end

end

# Serial implementation of SMC and cSMC
@inline function _smcSerial!(model::SMCModel, smcio::SMCIO{Particle},
  ref::Union{Nothing, Vector{Particle}} = nothing,
  refout::Union{Nothing, Vector{Particle}} = nothing) where Particle
  zetas = smcio.zetas
  zetaAncs = smcio.internal.zetaAncs
  lws = smcio.internal.lws
  ws = smcio.ws
  as = smcio.internal.as
  engine = getSMCRNG()
  pScratch = smcio.internal.particleScratch
  smcio.internal.nresamples = 0

  logZhats = smcio.logZhats
  lZ = 0.0

  _iota!(smcio.eves)

  for p = 1:smcio.n
    if p > 1
      _copyParticles!(zetaAncs, zetas, as)
      smcio.internal.oldEves .= smcio.eves
      _setEves!(smcio.eves, smcio.internal.oldEves, as)
    end
    if ref == nothing
      _mutateParticles!(zetas, engine, p, model.M!, zetaAncs, pScratch)
    else
      _mutateParticles!(zetas, engine, p, model.M!, zetaAncs, pScratch, ref[p])
    end

    _logWeightParticles!(lws, p, model.lG, zetas, pScratch,
      p == 1 || smcio.resample[p-1], ws)

    smcio.internal.maxlw = maximum(lws)
    ws .= exp.(lws .- smcio.internal.maxlw)
    smcio.internal.mws = mean(ws)
    lZ += smcio.internal.maxlw + log(smcio.internal.mws)
    ws ./= smcio.internal.mws
    @inbounds logZhats[p] = lZ
    @inbounds smcio.esses[p] = smcio.N / mapreduce(x -> x*x, +, ws)

    smcio.Vhat1s[p] = _Vhat1Serial(smcio)

    p < smcio.n && _resampleSerial!(smcio, p, ref != nothing)

    smcio.fullOutput && _intermediateOutputSerial!(smcio, p)
  end

  refout != nothing && SequentialMonteCarlo.pickParticle!(refout, smcio)
  return
end
