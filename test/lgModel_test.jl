using SequentialMonteCarlo
using SMCExamples.LinearGaussian
using Compat.Test

setSMCRNGs(0)

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
