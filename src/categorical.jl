import NonUniformRandomVariateGeneration.sampleSortedUniforms!

# sample nans categorical r.v.s where Pr(X=i+offset) ∝ p[i]
# required: length(Fs) = length(p), unifs[start+1:start+nans] are valid
@inline function sampleSerial!(vs::Vector{Int64}, rng::SMCRNG, p::V,
  nans::Int64, unifs::Vector{Float64}, Fs::V, start::Int64 = 0,
  offset::Int64 = 0) where V<:AbstractVector{Float64}
  # generateSortedUniformsBSerial!(unifs, start, nans, rng)
  sampleSortedUniforms!(unifs, start, nans, rng)
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
