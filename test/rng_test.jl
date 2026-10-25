using Random
using SequentialMonteCarlo

mutable struct RNGTestParticle
  x::Float64
  RNGTestParticle() = new(0.0)
end

@testset "Task-local RNG and migrating callbacks" begin
  function mutate!(dest, rng::AbstractRNG, p, ancestor, scratch)
    # A callback may yield and resume on another physical thread.
    yield()
    dest.x = (p == 1 ? 0.0 : ancestor.x) + randn(rng)
  end
  model = SMCModel(mutate!, (p, x, scratch) -> 0.0, 3,
    RNGTestParticle, Nothing)
  for chunks in unique([1, min(Threads.nthreads(), 4)])
    function runseeded()
      Random.seed!(123)
      io = SMCIO{RNGTestParticle, Nothing}(64, 3, chunks, true)
      smc!(model, io)
      return io
    end
    first = fetch(Threads.@spawn runseeded())
    second = fetch(Threads.@spawn runseeded())
    @test first.logZhats == zeros(3)
    @test first.esses == ones(3)
    @test [p.x for p in first.zetas] == [p.x for p in second.zetas]
    @test all(isfinite(p.x) for p in first.zetas)
    @test issorted(first.eves)
    ref = [RNGTestParticle() for _ in 1:3]
    for p in ref
      p.x = 9.0
    end
    fetch(Threads.@spawn csmc!(model, first, ref, ref))
    @test all(first.allZetas[p][1].x == 9.0 for p in 1:3)
    @test all(all(1 .<= a .<= 64) for a in first.allAs)
  end
end
