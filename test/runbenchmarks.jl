using SequentialMonteCarlo
import SMCExamples.LinearGaussian.defaultLGModel
import SMCExamples.MVLinearGaussian.defaultMVLGModel
import SMCExamples.SMCSampler.defaultSMCSampler
import SMCExamples.Lorenz96.defaultLorenzModel

using Dates
using BenchmarkTools, Compat

function benchMachineInfo()
  s::String = "$(Sys.MACHINE)\n" *
  "$(Sys.cpu_info()[1].model) ($(Sys.CPU_NAME)) ; " *
    "$(Sys.CPU_THREADS) Logical cores\n" *
    "Julia $VERSION using $(Threads.nthreads()) threads"
  return s
end


function dateandtime()
  datetimenow = now()
  s::String = Dates.format(datetimenow, "e, dd u yyyy HH:MM:SS")
  return s
end

function closestN(N::Int64, nt::Int64)
  mod(N, nt) == 0 && return N
  return round(Int64, N / nt) * nt
end

function bench(model, range)
  nt = Threads.nthreads()
  println("log₂N  Benchmark")
  for k in range
    print("$(lpad(string(k),2,string(0)))   ")
    N::Int64 = closestN(2^k, nt)
    smcio = SMCIO{model.particle, model.pScratch}(N, model.maxn, nt, false)
    @btime smc!($model, $smcio)
  end
end

println("Running Benchmarks: $(dateandtime())")
println(benchMachineInfo())

range = [10, 15, 20]

println("Linear Gaussian Model, n = 10")
model, theta, ys, ko = defaultLGModel(10)
bench(model, range)

println("Multivariate Linear Gaussian Model, d = 10, n = 10")
model, theta, ys, ko = defaultMVLGModel(10, 10)
bench(model, range)

println("SMC Sampler Example, n = 12")
model, ltarget = defaultSMCSampler()
bench(model, range)

println("Lorenz96 Example, d = 8, n = 10")
model, theta, ys = defaultLorenzModel(8, 10)
bench(model, range)

println("Finished: $(dateandtime())")
