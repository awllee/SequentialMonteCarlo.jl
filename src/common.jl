import NonUniformRandomVariateGeneration.sampleCategorical
import Compat: hasmethod, copyto!

# basic routines used by both smcSerial and smcParallel

# adapted from https://discourse.julialang.org/t/how-to-copy-all-fields-without-changing-the-referece/945/5
@inline @generated function particleCopy!(dest, src)
  fieldNames = fieldnames(dest)
  fieldTypes = dest.types
  numFields = length(fieldNames)
  expressions = Array{Expr}(undef, numFields)

  for i = 1:numFields
    fieldName = fieldNames[i]
    fieldType = fieldTypes[i]
    @assert !fieldType.mutable || hasmethod(copyto!, (fieldType,
      fieldType)) "$fieldName::$fieldType : copyto! must exist for mutable Particle fields"
    if fieldType.mutable
      @inbounds expressions[i] = :(copyto!(dest.$fieldName, src.$fieldName))
    else
      @inbounds expressions[i] = :(dest.$fieldName = src.$fieldName)
    end
  end
  body = Expr(:block, expressions...)

  quote
    # $(Expr(:meta, :inline))
    $body
    return
  end
end

## better not to inline ; enough work should be done in the loop
function _copyParticles!(out::Vec, in::Vec) where {Particle,
  Vec<:AbstractVector{Particle}}
  for i in eachindex(out)
    @inbounds particleCopy!(out[i], in[i])
  end
end

## better not to inline ; enough work should be done in the loop
function _copyParticles!(zetaAncs::Vec1, zetas::Vector{Particle},
  as::Vec2) where {Particle, Vec1<:AbstractVector{Particle},
    Vec2<:AbstractVector{Int64}}
  for i in eachindex(zetaAncs)
    @inbounds particleCopy!(zetaAncs[i], zetas[as[i]])
  end
end

## better not to inline ; enough work should be done in the loop
function _setEves!(eves::Vec1, oldEves::Vec2, as::Vec1) where
  {Vec1<:AbstractVector{Int64}, Vec2<:AbstractVector{Int64}}
  for i in eachindex(eves)
    @inbounds eves[i] = oldEves[as[i]]
  end
end

## better not to inline ; enough work should be done in the loop
function _mutateParticles!(zetas::Vec, engine::RNG, p::Int64, M!::F,
  zetaAncs::Vec, pScratch::ParticleScratch) where {Particle,
  Vec<:AbstractVector{Particle}, F<:Function, ParticleScratch}
  for j in eachindex(zetas)
    @inbounds M!(zetas[j], engine, p, zetaAncs[j], pScratch)
  end
end

## better not to inline ; enough work should be done in the loop
function _mutateParticles!(zetas::Vec, engine::RNG, p::Int64, M!::F,
  zetaAncs::Vec, pScratch::ParticleScratch, xref::Particle) where {Particle,
  Vec<:AbstractVector{Particle}, F<:Function, ParticleScratch}
  particleCopy!(zetas[1], xref)
  for j in 2:length(zetas)
    @inbounds M!(zetas[j], engine, p, zetaAncs[j], pScratch)
  end
end

## better not to inline ; enough work should be done in the loop
function _logWeightParticles!(lws::Vec1, p::Int64, lG::F,
  zetas::Vec2, pScratch::ParticleScratch, resample::Bool, oldWs::Vec1) where
  {Vec1<:AbstractVector{Float64}, F<:Function, Particle,
  Vec2<:AbstractVector{Particle}, ParticleScratch}
  for j in eachindex(zetas)
    @inbounds lws[j] = lG(p, zetas[j], pScratch)
  end
  !resample && (lws .+= log.(oldWs))
end

@inline function pickParticle!(path::Vector{Particle}, smcio::SMCIO{Particle}) where Particle
  @assert smcio.fullOutput
  @assert length(path) == smcio.n
  n::Int64 = smcio.n
  allAs::Vector{Vector{Int64}} = smcio.allAs
  allZetas::Vector{Vector{Particle}} = smcio.allZetas
  rng::RNG = getRNG()
  k::Int64 = sampleCategorical(smcio.ws, rng)

  @inbounds particleCopy!(path[n], allZetas[n][k])
  for p = n-1:-1:1
    @inbounds k = allAs[p][k]
    @inbounds particleCopy!(path[p], allZetas[p][k])
  end
  return
end

@inline function pickParticleBS!(path::Vector{Particle},
  smcio::SMCIO{Particle, ParticleScratch}, lM::F) where {Particle,
  ParticleScratch, F<:Function}
  @assert smcio.fullOutput
  @assert length(path) == smcio.n
  n::Int64 = smcio.n
  N::Int64 = smcio.N
  allAs::Vector{Vector{Int64}} = smcio.allAs
  allZetas::Vector{Vector{Particle}} = smcio.allZetas
  allWs::Vector{Vector{Float64}} = smcio.allWs
  bws::Vector{Float64} = smcio.internal.scratch1
  pScratch::ParticleScratch = smcio.internal.particleScratch

  rng = getRNG()
  k::Int64 = sampleCategorical(smcio.ws, rng)

  @inbounds particleCopy!(path[n], allZetas[n][k])
  for p = n-1:-1:1
    @inbounds bws .= log.(allWs[p])
    for j in 1:N
      @inbounds bws[j] += lM(p+1, allZetas[p][j], path[p+1], pScratch)
    end
    m::Float64 = maximum(bws)
    bws .= exp.(bws .- m)
    k = sampleCategorical(bws, rng)
    @inbounds particleCopy!(path[p], allZetas[p][k])
  end
end
