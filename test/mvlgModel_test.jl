using SequentialMonteCarlo
using RNGPool
using SMCExamples.MVLinearGaussian
using StaticArrays
import SMCExamples.Particles.MVFloat64Particle
using Compat.Test

setRNGs(0)

function mvlgtest(d::Int64, nthreads::Int64)
  model, theta, ys, ko = MVLinearGaussian.defaultMVLGModel(5, d)

  numParticles = 2^16

  smcio = SMCIO{model.particle, model.pScratch}(numParticles, model.maxn,
    nthreads, false)
  smc!(model, smcio)

  @test smcio.logZhats ≈ ko.logZhats atol=0.1*d
end

function testmvlgcsmc(nthreads::Int64, essThreshold::Float64)
  model, theta, ys, ko = MVLinearGaussian.defaultMVLGModel(2, 10)
  nsamples = 2^12

  smcio = SMCIO{model.particle, model.pScratch}(256, model.maxn, nthreads, true,
    essThreshold)

  v = Vector{MVFloat64Particle{2}}(undef, 10)
  for p = 1:10
    v[p] = MVFloat64Particle{2}()
    v[p].x .= zeros(MVector{2, Float64})
  end

  meanEstimates = Vector{MVector{2,Float64}}(undef, 10)
  for p = 1:10
    meanEstimates[p] = zeros(MVector{2, Float64})
  end
  for i = 1:nsamples
    csmc!(model, smcio, v, v)
    for p = 1:10
      meanEstimates[p] .+= v[p].x
    end
  end
  meanEstimates ./= nsamples

  testapproxequal((x -> x[1]).(meanEstimates), (x -> x[1]).(ko.smoothingMeans),
    0.1, false)
  testapproxequal((x -> x[2]).(meanEstimates), (x -> x[2]).(ko.smoothingMeans),
    0.1, false)
end

mvlgtest(1, 1)
mvlgtest(1, Threads.nthreads())
mvlgtest(2, 1)
mvlgtest(2, Threads.nthreads())

testmvlgcsmc(1, 2.0)
testmvlgcsmc(1, 0.5)
nt = min(Threads.nthreads(), 4)
testmvlgcsmc(nt, 2.0)
testmvlgcsmc(nt, 0.5)
