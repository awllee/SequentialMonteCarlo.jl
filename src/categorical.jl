# sample Binomial r.v. using inverse transform
@inline function sampleInverseCDFBinom(n::Int64, p::Float64, rng::SMCRNG)
  q::Float64 = 1-p
  s::Float64 = p/q
  k::Int64 = 0
  f::Float64 = q^n
  u::Float64 = rand(rng)
  while u > f
    u -= f
    k += 1
    f *= (n + 1 - k) / k * s
  end
  return k
end

let
  btrd_fc_table = [
    0.08106146679532726,
    0.04134069595540929,
    0.02767792568499834,
    0.02079067210376509,
    0.01664469118982119,
    0.01387612882307075,
    0.01189670994589177,
    0.01041126526197209,
    0.009255462182712733,
    0.008330563433362871]
  global @inline function btrd_fc(k::Int64)
    if k < 10
      return btrd_fc_table[k+1]
    else
      tmp::Float64 = 1/(k+1)
      return (1/12 - (1/360 - tmp * tmp / 1260) * tmp * tmp) * tmp
    end
  end
end

## Generate a Binomial r.v. using the algorithm of
## W. Hörmann. The Generation of Binomial Random Variates, Journal of
## Statistical Computation and Simulation 46, 1993.
## intended for use when n*p >= 10 and p <= 0.5
@inline function btrd(n::Int64, p::Float64, rng::SMCRNG)
  # 0. setup
  m::Int64 = floor(Int64, (n + 1) * p)
  r::Float64 = p / (1 - p)
  nr::Float64 = (n + 1) * r
  npq::Float64 = n * p * (1 - p)
  sqrt_npq::Float64 = sqrt(npq)
  b::Float64 = 1.15 + 2.53 * sqrt_npq
  a::Float64 = -0.0873 + 0.0248 * b + 0.01 * p
  c::Float64 = n * p + 0.5
  α::Float64 = (2.83 + 5.1 / b) * sqrt_npq
  vr::Float64 = 0.92 - 4.2 / b
  ur_vr::Float64 = 0.86 * vr
  u::Float64 = 0.0
  v::Float64 = 0.0
  us::Float64 = 0.0
  k::Int64 = 0
  km::Int64 = 0
  f::Float64 = 0.0
  i::Int64 = 0
  while true
    # 1.
    v = rand(rng)
    if v <= ur_vr
      u = v / vr - 0.43
      return floor(Int64, (2 * a / (0.5 - abs(u)) + b) * u + c)
    end
    # 2.
    if v >= vr
      u = rand(rng) - 0.5
    else
      u = v / vr - 0.93
      u = sign(u) * 0.5 - u
      v = vr * rand(rng)
    end
    # 3.0
    us = 0.5 - abs(u)
    k = floor(Int64, (2 * a / us + b) * u + c)
    if k < 0 || k > n
      continue
    end
    v = v * α / (a / (us * us) + b)
    km = abs(k - m)
    # [ if km > 15.0 goto 3.2 ]
    if km <= 15.0
      # 3.1 [ recursive evaluation of f(k) ]
      f = 1.0
      i = min(m, k)
      if m < k
        while i != k
          i += 1
          f *= (nr / i - r)
        end
      elseif m > k
        while i != m
          i += 1
          v *= (nr / i - r)
        end
      end
      if v <= f
        return k
      else
        continue
      end
    else
      # 3.2
      v = log(v)
      ρ = (km / npq) * (((km / 3 + 0.625) * km + 1/6) / npq + 0.5)
      t = - km * km / (2 * npq)
      if v < t - ρ
        return k
      end
      if v > t + ρ
        continue
      end
      # 3.3
      nm = n - m + 1
      h = (m + 0.5) * log((m + 1) / (r * nm)) + btrd_fc(m) + btrd_fc(n - m)
      # 3.4
      nk = n - k + 1
      if v <= h + (n+1)*log(nm/nk)+(k+0.5)*log(nk*r/(k+1))-btrd_fc(k)-btrd_fc(n-k)
        return k
      end
    end
  end
end

# choose which algorithm to use for given (n,p)
@inline function sampleBinomial(n::Int64, p::Float64, rng::SMCRNG)
  if n * p >= 10 && n * (1.0 - p) >= 10
    if p <= 0.5
      return btrd(n, p, rng)
    else
      return n - btrd(n, 1.0 - p, rng)
    end
  else
    if p <= 0.5
      return sampleInverseCDFBinom(n, p, rng)
    else
      return n - sampleInverseCDFBinom(n, 1.0 - p, rng)
    end
  end
end

# sample a multinomial(n,p) r.v. using binomial r.v.s
@inline function sampleMultinomial!(n::Int64, p::Vector{Float64},
  result::Vector{Int64}, rng::SMCRNG)
  nleft::Int64 = n
  tmp::Float64 = 0.0
  for i = 1:length(p)-1
    @inbounds result[i] = sampleBinomial(nleft, min(1.0, p[i] / (1.0 - tmp)), rng)
    @inbounds nleft -= result[i]
    @inbounds tmp += p[i]
  end
  @inbounds result[length(p)] = nleft;
  return
end

# generate a sorted array of uniform(0,1) r.v.s
# This is the uniform spacing method, Algorithm B on p. 214 of
# L. Devroye. Non-uniform random variate generation. 1986.
@inline function generateSortedUniformsBSerial!(vec::Vector{Float64},
  start::Int64, n::Int64, rng::SMCRNG)
  # randexp!(rng, vec)
  s::Float64 = 0.0
  for i = 1:n
    @inbounds vec[start + i] = randexp(rng)
    @inbounds s += vec[start + i]
  end
  G::Float64 = 1.0 / (s + randexp(rng))
  @inbounds vec[start + 1] *= G
  for i = 2:n
    @inbounds vec[start + i] *= G
    @inbounds vec[start + i] += vec[start + i - 1]
  end
end

# sample nans categorical r.v.s where Pr(X=i+offset) ∝ p[i]
# required: length(Fs) = length(p), unifs[start+1:start+nans] are valid
@inline function sampleSerial!(vs::Vector{Int64}, rng::SMCRNG, p::V,
  nans::Int64, unifs::Vector{Float64}, Fs::V, start::Int64 = 0,
  offset::Int64 = 0) where V<:AbstractVector{Float64}
  generateSortedUniformsBSerial!(unifs, start, nans, rng)
  cumsum!(Fs, p)
  @inbounds maxFs = Fs[length(p)]
  j = 1
  for i = 1:nans
    @inbounds while maxFs * unifs[start + i] > Fs[j]
      j += 1
    end
    @inbounds vs[start + i] = j + offset
  end
end

@inline function sampleOne(ws::Vector{Float64}, rng::SMCRNG)
  s::Float64 = sum(ws)
  u::Float64 = rand(rng) * s
  i::Int64 = 1
  v::Float64 = ws[1]
  while u > v
    i += 1
    @inbounds v += ws[i]
  end
  return i
end
