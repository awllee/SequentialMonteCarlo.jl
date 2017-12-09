using SequentialMonteCarlo
import SMCExamples.SMCSampler.defaultSMCSampler
using Compat.Test

setSMCRNGs(0)

model, ltarget = defaultSMCSampler()

numParticles = 1024*256

smcio = SMCIO{model.particle, model.pScratch}(numParticles, model.maxn, 1,
  false)
smc!(model, smcio)

@test smcio.logZhats[model.maxn] ≈ 0.0 atol=0.1

smcio = SMCIO{model.particle, model.pScratch}(numParticles, model.maxn,
  Threads.nthreads(), false)
smc!(model, smcio)

@test smcio.logZhats[model.maxn] ≈ 0.0 atol=0.1
