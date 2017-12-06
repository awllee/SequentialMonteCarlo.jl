using SequentialMonteCarlo
import SMCExamples.LinearGaussian.defaultLGModel
import SMCExamples.MVLinearGaussian.defaultMVLGModel
import SMCExamples.SMCSampler.defaultSMCSampler
import SMCExamples.Lorenz96.defaultLorenzModel

VERSION.minor > 6 && using Dates

using BenchmarkTools, Compat
import Hwloc

function benchMachineInfo()
  s::String = "$(Sys.MACHINE) ; " *
    "$(Sys.cpu_info()[1].model) ($(Sys.cpu_name)) ; " *
    "$(Hwloc.num_physical_cores()) Physical, $(Sys.CPU_CORES) Logical"
  return s
end

function dateandtime()
  datetimenow = now()
  s::String = Dates.format(datetimenow, "e, dd u yyyy HH:MM:SS")
  return s
end

function bench(model, range)
  nt = Threads.nthreads()
  println("log₂N  Threads  Benchmark")
  for k in range
    print("$(lpad(k,2,0))     $(lpad(nt,2,0))     ")
    smcio = SMCIO{model.particle, model.pScratch}(2^k, model.maxn, nt, false)
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
bench(model, range[1:2])

println("Finished: $(dateandtime())")
