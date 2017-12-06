using Documenter, SequentialMonteCarlo

makedocs(
  modules = [SequentialMonteCarlo],
  format = :html,
  sitename = "SequentialMonteCarlo.jl",
  authors = "Anthony Lee",
  pages = [
    "Contents" => "contents.md"
    "Introduction" => "intro.md"
    "SMC integrals" => "smcintegrals.md"
    "SMC algorithm" => "smcalgorithm.md"
    "Theoretical properties" => "smctheory.md"
    "Variance estimators" => "smcve.md"
    "Adaptive resampling" => "smcadaptive.md"
    "Conditional SMC" => "csmc.md"
    "Implementation notes" => "impl.md"
    "SMC interface" => "smcinterface.md"
    "Types and functions" => "guide.md"
    "References" => "refs.md"
  ]
)
