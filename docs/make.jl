using Documenter, SequentialMonteCarlo

makedocs(
  modules = [SequentialMonteCarlo],
  format = :html,
  sitename = "SequentialMonteCarlo.jl",
  authors = "Anthony Lee",
  pages = [
    "Contents" => "index.md"
    "Introduction" => "intro.md"
    "SMC integrals" => "smcintegrals.md"
    "SMC algorithm" => "smcalgorithm.md"
    "Theoretical properties" => "smctheory.md"
    "Variance estimators" => "smcve.md"
    "Adaptive resampling" => "smcadaptive.md"
    "Conditional SMC" => "csmc.md"
    "Implementation notes" => "impl.md"
    "Performance tips" => "performance.md"
    "SMC interface" => "smcinterface.md"
    "Types and functions" => "guide.md"
    "Hidden Markov models" => "hmm.md"
    "Benchmarks" => "bench.md"
    "References" => "refs.md"
  ]
)

deploydocs(
  repo = "github.com/awllee/SequentialMonteCarlo.jl.git",
  target = "build",
  julia  = "1.0",
  osname = "linux",
  deps = nothing,
  make = nothing,
)
