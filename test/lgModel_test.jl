using SequentialMonteCarlo
using RNGPool
using SMCExamples.LinearGaussian
using Test

setRNGs(0)

model, theta, ys, ko = LinearGaussian.defaultLGModel(5)

numParticles = 2^16

smcio = SMCIO{model.particle, model.pScratch}(numParticles, model.maxn, 1,
  false)
smc!(model, smcio)

@test smcio.logZhats ≈ ko.logZhats atol=0.1

nthreads = Threads.nthreads()
if nthreads > 1
  smcio = SMCIO{model.particle, model.pScratch}(numParticles, model.maxn,
    nthreads, false)
  smc!(model, smcio)

  @test smcio.logZhats ≈ ko.logZhats atol=0.1
end

lM = LinearGaussian.makelM(theta)
smcio = SMCIO{model.particle, model.pScratch}(numParticles, model.maxn,
  nthreads, true)
smc!(model, smcio)
path = Vector{model.particle}(undef, smcio.n)
for i in 1:smcio.n
  path[i] = model.particle()
end
SequentialMonteCarlo.pickParticleBS!(path, smcio, lM)
