import SequentialMonteCarlo.sampleMultinomial!
import SequentialMonteCarlo.setSMCRNGs
import SequentialMonteCarlo.getSMCRNG
using Compat.Test

setSMCRNGs(0)

function multinomial_test()
  p = [.2, .5, .3]
  n = 32
  rng = getSMCRNG()

  m = 2^18
  array = zeros(Float64, 3)
  array2 = zeros(Float64, 3)
  result = Vector{Int64}(uninitialized, 3)
  E = n*p
  V = n*p.*(1.0 .- p)

  for i=1:m
    sampleMultinomial!(n, p, result, rng)
    array .+= result - E
    array2 .+= (result - E).^2
  end

  @test maximum(abs.(array / m - [0.0, 0.0, 0.0])) < 0.01
  @test maximum(abs.(array2 / m - V)) < 0.05
end

multinomial_test()
