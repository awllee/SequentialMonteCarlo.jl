"""
    smc!(model::SMCModel, smcio::SMCIO)
Run the SMC algorithm for the given model and input/output arguments.

If ```smcio.nthreads == 1``` the algorithm will run in serial.
"""
function smc!(model::SMCModel, smcio::SMCIO)
  if smcio.nthreads > 1
    _smcParallel!(model, smcio)
  else
    _smcSerial!(model, smcio)
  end
end

"""
    csmc!(model::SMCModel, smcio::SMCIO, ref::Vector{Particle}, refout::Vector{Particle})

Run the conditional SMC algorithm for the given model, input/output arguments,
reference path and output path.

It is permitted for ref and refout to be the same. If ```smcio.nthreads == 1```
the algorithm will run in serial.
"""
function csmc!(model::SMCModel, smcio::SMCIO{Particle}, ref::Vector{Particle},
  refout::Vector{Particle}) where Particle
  if smcio.nthreads > 1
    _smcParallel!(model, smcio, ref, refout)
  else
    _smcSerial!(model, smcio, ref, refout)
  end
end

@inline function _getZetas(smcio::SMCIO{Particle}, hat::Bool, p::Int64) where
  Particle
  @assert 1 <= p <= smcio.n
  p == smcio.n && return smcio.zetas
  @assert smcio.fullOutput "full SMC output required for this operation"
  return smcio.allZetas[p]
end

@inline function _getWs(smcio::SMCIO, hat::Bool, p::Int64)
  @assert 1 <= p <= smcio.n
  if hat
    p == smcio.n && return smcio.ws
    @assert smcio.fullOutput "full SMC output required for this operation"
    return smcio.allWs[p]
  else
    if p == 1 || smcio.resample[p-1] == true
      return smcio.internal.vecOnes
    else
      @assert smcio.fullOutput "full SMC output required for this operation"
      return smcio.allWs[p-1]
    end
  end
end

"""
    eta(smcio::SMCIO{Particle}, f::F, hat::Bool, p::Int64) where {Particle, F<:Function}

Compute:
  - ```!hat```: ``\\eta^N_p(f)``
  - ```hat```:  ``\\hat{\\eta}_p^N(f)``
"""
function eta(smcio::SMCIO{Particle}, f::F, hat::Bool, p::Int64) where
  {Particle, F<:Function}
  zetas::Vector{Particle} = _getZetas(smcio, hat, p)
  ws::Vector{Float64} = _getWs(smcio, hat, p)
  v = f(zetas[1]) * ws[1]
  for i = 2:smcio.N
    @inbounds v += f(zetas[i]) * ws[i]
  end
  return v / smcio.N
end

"""
    allEtas(smcio::SMCIO, f::F, hat::Bool) where F<:Function

Compute ```eta(smcio::SMCIO, f::F, hat::Bool, p)``` for p in {1, …, smcio.n}
"""
function allEtas(smcio::SMCIO, f::F, hat::Bool) where F<:Function
  T = typeof(f(smcio.zetas[1]) / smcio.N)
  result::Vector{T} = Vector{T}(undef, smcio.n)
  for p = 1:smcio.n
    @inbounds result[p] = eta(smcio, f, hat, p)
  end
  return result
end

"""
    slgamma(smcio::SMCIO, f::F, hat::Bool, p::Int64) where {Particle, F<:Function}

Compute:
- ```!hat```: ``(\\eta^N_p(f) \\geq 0, \\log |\\gamma^N_p(f)|)``
- ```hat```:  ``(\\hat{\\eta}^N_p(f) \\geq 0, \\log |\\hat{\\gamma}_p^N(f)|)``
The result is returned as a ```Tuple{Bool, Float64}```: the first component
represents whether the returned value is non-negative, the second is the
logarithm of the absolute value of the approximation.
"""
function slgamma(smcio::SMCIO, f::F, hat::Bool, p::Int64) where
  F<:Function
  @assert 1 <= p <= smcio.n "p must be between 1 and smcio.n"
  idx::Int64 = p - 1 + hat
  logval = idx == 0 ? 0.0 : smcio.logZhats[idx]
  v::Float64 = eta(smcio, f, hat, p)
  logval += log(abs(v))
  return (v >= 0, logval)
end

"""
    allGammas(smcio::SMCIO, f::F, hat::Bool) where F<:Function

Compute ```slgamma(smcio::SMCIO, f::F, hat::Bool, p)``` for p in {1, …, smcio.n}
"""
function allGammas(smcio::SMCIO, f::F, hat::Bool) where F<:Function
  result::Vector{Tuple{Bool, Float64}} =
    Vector{Tuple{Bool, Float64}}(undef, smcio.n)
  for p = 1:smcio.n
    @inbounds result[p] = slgamma(smcio, f, hat, p)
  end
  return result
end

function _centreFunction(smcio::SMCIO, f::F, hat::Bool,
  p::Int64) where F<:Function
  mu = eta(smcio, f, hat, p)
  function cf(x)
    return f(x) - mu
  end
  return cf
end

"""
    V(smcio::SMCIO{Particle}, f::F, hat::Bool, centred::Bool, p::Int64) where
      {Particle, F<:Function}
Compute:
- ```!hat & !centred```: ``V^N_p(f)``
- ```!hat & centred```: ``V^N_p(f-\\eta_p^N(f))``
- ```hat & !centred```: ``\\hat{V}_p^N(f)``
- ```hat & centred```:  ``\\hat{V}_p^N(f-\\hat{\\eta}_p^N(f))``
"""
function V(smcio::SMCIO{Particle}, f::F, hat::Bool, centred::Bool,
  p::Int64) where {Particle, F<:Function}

  centred && return V(smcio, _centreFunction(smcio, f, hat, p), hat, false, p)

  zetas::Vector{Particle} = _getZetas(smcio, hat, p)
  ws::Vector{Float64} = _getWs(smcio, hat, p)
  if p == smcio.n
    eves = smcio.eves
  else
    eves = smcio.allEves[p]
  end

  etafSq = eta(smcio, f, hat, p)^2
  N::Int64 = smcio.N
  i::Int64 = 1
  v::Float64 = 0.0
  while i <= N
    @inbounds e = eves[i]
    tmp::Float64 = 0.0
    @inbounds while i <= N && eves[i] == e
      @inbounds tmp += ws[i] * f(zetas[i])
      i += 1
    end
    tmp /= N
    v += tmp*tmp
  end
  q::Int64 = 1
  for i = 1:p-1
    q += smcio.resample[i]
  end
  return etafSq - (etafSq-v)*(N/(N-1))^q
end

## currently single-threaded only
"""
    vpns(smcio::SMCIO, f::F, hat::Bool, centred::Bool, n::Int64) where F<:Function

Compute a vector of the values of, for p = 1, …, n,
- ```!hat & !centred```: ``v^N_{p,n}(f)``
- ```!hat & centred```:  ``v^N_{p,n}(f-\\eta_n^N(f))``
- ```hat & !centred```:  ``\\hat{v}_{p,n}^N(f)``
- ```hat & centred```:   ``\\hat{v}_{p,n}^N(f-\\hat{\\eta}_n^N(f))``

Note: if ```essThreshold <= 1.0```, and resampling did not occur at every time,
the length of the output will be less than n.
"""
function vpns(smcio::SMCIO, f::F, hat::Bool, centred::Bool,
  n::Int64) where F<:Function

  centred && return vpns(smcio, _centreFunction(smcio, f, hat, n), hat, false, n)

  N = smcio.N
  as = smcio.allAs
  eves = smcio.allEves
  ws = smcio.allWs
  zetas = smcio.allZetas[n]

  etafSq = eta(smcio, f, hat, n)^2

  Slocal = smcio.internal.scratch1
  if hat
    for i = 1:N
      @inbounds Slocal[i] = f(zetas[i]) * ws[n][i] / N
    end
  else
    if n > 1 && !smcio.resample[n-1]
      for i = 1:N
        @inbounds Slocal[i] = f(zetas[i]) * ws[n-1][i] / N
      end
    else
      for i = 1:N
        @inbounds Slocal[i] = f(zetas[i]) / N
      end
    end
  end
  Stmp = smcio.internal.scratch2
  Glocal = smcio.internal.scratch2
  fill!(Glocal, 0.0)
  if n > 1
    for i = 1:N
      @inbounds Glocal[eves[n-1][i]] += ws[n-1][i] / N
    end
  end

  result = Vector{Float64}(undef, n)
  factor0::Float64 = (N/(N-1))^n
  factor1::Float64 = (N-1)*factor0
  vp::Float64 = 0.0
  for i = 1:N
    @inbounds vp += Slocal[i]*Slocal[i] * (1.0 - Glocal[eves[n][i]])
  end
  result[n] = factor1 * vp

  for p = n-1:-1:1
    # Slocal = S[p+1,:]
    # Glocal = G[p-1,:]
    fill!(Glocal, 0.0)
    if p != 1
      for i = 1:N
        @inbounds Glocal[eves[p-1][i]] += ws[p-1][i] / N
      end
    end
    i::Int64 = 1
    vp = 0.0
    while i <= N
      @inbounds a::Int64 = as[p][i]
      R1::Float64 = 0.0
      R2::Float64 = 0.0
      @inbounds while i <= N && as[p][i] == a
        @inbounds R1 += Slocal[i]
        @inbounds R2 += Slocal[i]*Slocal[i]
        i += 1
      end
      @inbounds vp += (R1*R1 - R2) * (1.0 - Glocal[eves[p][a]])
    end
    @inbounds result[p] = factor1 * vp
    copyto!(Stmp, Slocal)
    fill!(Slocal, 0.0)
    for i = 1:N
      @inbounds Slocal[as[p][i]] += Stmp[i]
    end
  end

  vp = 0.0
  for i = 1:N
    @inbounds vp += Slocal[i]*Slocal[i]
  end
  mstar::Float64 = factor0 * (etafSq - vp)
  result .-= mstar

  actualResult = Vector{Float64}(undef, sum(smcio.resample[1:n-1]) + 1)
  j::Int64 = 0
  for i = 1:length(actualResult)-1
    j += 1
    while !smcio.resample[j]
      j += 1
    end
    actualResult[i] = result[j]
  end
  actualResult[length(actualResult)] = result[n]
  return actualResult
end

"""
    v(smcio::SMCIO, f::F, hat::Bool, centred::Bool, n::Int64) where F<:Function

Compute:
- ```!hat & !centred```: ``v^N_n(f)``
- ```!hat & centred```: ``v^N_n(f-\\eta_n^N(f))``
- ```hat & !centred```: ``\\hat{v}_n^N(f)``
- ```hat & centred```:  ``\\hat{v}_n^N(f-\\hat{\\eta}_n^N(f))``
"""
function v(smcio::SMCIO, f::F, hat::Bool, centred::Bool,
  n::Int64) where F<:Function
  return sum(vpns(smcio, f, hat, centred, n))
end
