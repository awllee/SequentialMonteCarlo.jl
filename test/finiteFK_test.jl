using SequentialMonteCarlo
using SMCExamples.FiniteFeynmanKac
import SMCExamples.Particles.Int64Particle
using Compat.Test

function Id(p::Int64Particle)
  return p.x
end

function One(p::Int64Particle)
  return 1.0
end

function IdMinus2(p::Int64Particle)
  return p.x - 2
end

function testEtas(smcio, ffkout, verbose::Bool)
  testapproxequal(
    SequentialMonteCarlo.allEtas(smcio, Id, false),
    FiniteFeynmanKac.allEtas(ffkout, Id, false), 0.01, verbose)
  testapproxequal(
    SequentialMonteCarlo.allEtas(smcio, Id, true),
    FiniteFeynmanKac.allEtas(ffkout, Id, true), 0.01, verbose)
end

function testZs(smcio, ffkout, verbose::Bool)
  testapproxequal(smcio.logZhats, ffkout.logZhats, 0.1, verbose)
end


function testapproxequal(x::Vector{Tuple{Bool,Float64}},
  y::Vector{Tuple{Bool,Float64}}, tol::Float64, verbose::Bool)
  xNonneg = (p->p[1]).(x)
  xLogVals = (p->p[2]).(x)
  yNonneg = (p->p[1]).(y)
  yLogVals = (p->p[2]).(y)
  @test xNonneg == yNonneg
  testapproxequal(xLogVals, yLogVals, tol, verbose)
end

function testGammas(smcio, ffkout, verbose::Bool)
  testapproxequal(SequentialMonteCarlo.allGammas(smcio, IdMinus2, false),
    FiniteFeynmanKac.allGammas(ffkout, IdMinus2, false), 0.1, verbose)
  testapproxequal(SequentialMonteCarlo.allGammas(smcio, IdMinus2, true),
    FiniteFeynmanKac.allGammas(ffkout, IdMinus2, true), 0.1, verbose)
end

function testVhat1s(smcio, ffkout, ffk, verbose::Bool)
  testapproxequal(smcio.N * smcio.Vhat1s,
    FiniteFeynmanKac.allavarhat1s(ffk, ffkout, smcio.resample), 0.1, verbose)
end

function _testV(smcio, ffkout, ffk, f, p, tol, verbose::Bool)
  resample = smcio.resample
  testapproxequal(smcio.N * SequentialMonteCarlo.V(smcio, f, false, false, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, false, false, resample, p), tol,
    verbose)

  testapproxequal(smcio.N * SequentialMonteCarlo.V(smcio, f, false, true, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, false, true, resample, p), tol,
    verbose)

  testapproxequal(smcio.N * SequentialMonteCarlo.V(smcio, f, true, false, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, true, false, resample, p), tol,
    verbose)

  testapproxequal(smcio.N * SequentialMonteCarlo.V(smcio, f, true, true, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, true, true, resample, p), tol,
    verbose)
end

function testV(smcio, ffkout, ffk, verbose::Bool)
  # verbose = true
  tol::Float64 = 0.25
  _testV(smcio, ffkout, ffk, One, smcio.n, tol, verbose)
  for p = 1:smcio.n
    _testV(smcio, ffkout, ffk, Id, p, tol, verbose)
  end
end

function _testv(smcio, ffkout, ffk, f, p, tol, verbose::Bool)
  resample = smcio.resample
  testapproxequal(SequentialMonteCarlo.v(smcio, f, false, false, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, false, false, resample, p), tol,
    verbose)
  testapproxequal(SequentialMonteCarlo.v(smcio, f, false, true, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, false, true, resample, p), tol,
    verbose)
  testapproxequal(SequentialMonteCarlo.v(smcio, f, true, false, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, true, false, resample, p), tol,
    verbose)
  testapproxequal(SequentialMonteCarlo.v(smcio, f, true, true, p),
    FiniteFeynmanKac.avar(ffk, ffkout, f, true, true, resample, p), tol,
    verbose)
end

function testv(smcio, ffkout, ffk, verbose::Bool)
  tol::Float64 = 0.25
  _testv(smcio, ffkout, ffk, One, smcio.n, tol, verbose)
  for p = 1:smcio.n
    _testv(smcio, ffkout, ffk, Id, p, tol, verbose)
  end
end

function _testvpn(smcio, ffkout, ffk, f, p, tol, verbose::Bool)
  resample = smcio.resample
  testapproxequal(SequentialMonteCarlo.vpns(smcio, f, false, false, p),
    FiniteFeynmanKac.vpns(ffk, ffkout, f, false, false, resample, p), tol,
    verbose)
  testapproxequal(SequentialMonteCarlo.vpns(smcio, f, false, true, p),
    FiniteFeynmanKac.vpns(ffk, ffkout, f, false, true, resample, p), tol,
    verbose)
  testapproxequal(SequentialMonteCarlo.vpns(smcio, f, true, false, p),
    FiniteFeynmanKac.vpns(ffk, ffkout, f, true, false, resample, p), tol,
    verbose)
  testapproxequal(SequentialMonteCarlo.vpns(smcio, f, true, true, p),
    FiniteFeynmanKac.vpns(ffk, ffkout, f, true, true, resample, p), tol,
    verbose)
end

function testvpn(smcio, ffkout, ffk, verbose::Bool)
  tol::Float64 = 0.2
  _testvpn(smcio, ffkout, ffk, One, smcio.n, tol, verbose)
  for p = 1:smcio.n
    _testvpn(smcio, ffkout, ffk, Id, p, tol, verbose)
  end
end

function testFullOutput(ffk::FiniteFeynmanKac.FiniteFK, nthreads::Int64,
  essThreshold::Float64, verbose::Bool)
  ffkout = FiniteFeynmanKac.calculateEtasZs(ffk)
  model = FiniteFeynmanKac.makeSMCModel(ffk)
  n = model.maxn
  smcio = SMCIO{model.particle, model.pScratch}(2^20, n, nthreads, true,
    essThreshold)
  smc!(model, smcio)

  @testset "Finite FK: eta" begin
    testEtas(smcio, ffkout, verbose)
  end

  @testset "Finite FK: Zhat" begin
    testZs(smcio, ffkout, verbose)
  end

  @testset "Finite FK: slgamma" begin
    testGammas(smcio, ffkout, verbose)
  end

  @testset "Finite FK: Vhat1" begin
    testVhat1s(smcio, ffkout, ffk, verbose)
  end

  @testset "Finite FK: V" begin
    testV(smcio, ffkout, ffk, verbose)
  end

  @testset "Finite FK: vpn" begin
    testvpn(smcio, ffkout, ffk, verbose)
  end

  @testset "Finite FK: v" begin
    testv(smcio, ffkout, ffk, verbose)
  end
end

function _getFreqs(model::SMCModel, smcio::SMCIO, m::Int64, d::Int64)
  v::Vector{Int64Particle} = FiniteFeynmanKac.Int642Path(1, d, smcio.n)
  counts = zeros(Int64, d^smcio.n)
  result = Vector{Float64}(d^smcio.n)
  for i = 1:m
    csmc!(model, smcio, v, v)
    counts[FiniteFeynmanKac.Path2Int64(v, d)] += 1
  end
  result .= counts ./ m
  return result
end

function testcsmc(nthreads::Int64, essThreshold::Float64)
  d = 3
  n = 3
  ffk = FiniteFeynmanKac.randomFiniteFK(d, n)

  model = FiniteFeynmanKac.makeSMCModel(ffk)
  densities = Vector{Float64}(d^n)
  for i = 1:length(densities)
    densities[i] = FiniteFeynmanKac.fullDensity(ffk, FiniteFeynmanKac.Int642Path(i, d, n))
  end
  densities ./= sum(densities)

  nsamples = 2^12

  smcio = SMCIO{model.particle, model.pScratch}(max(4,nthreads), model.maxn,
    nthreads, true, essThreshold)
  freqs = _getFreqs(model, smcio, nsamples, d)

  testapproxequal(freqs, densities, 0.05, false)
end

setSMCRNGs(0)

verbose = false

d = 3
ffk = FiniteFeynmanKac.randomFiniteFK(d, 10)

println("\nSerial: always resampling\n")
@time testFullOutput(ffk, 1, 2.0, verbose)
println("\nSerial: adaptive resampling, essThreshold = 0.75\n")
@time testFullOutput(ffk, 1, 0.75, verbose)
println("\nSerial: adaptive resampling, essThreshold = 0.25\n")
@time testFullOutput(ffk, 1, 0.25, verbose)

println("\nParallel: always resampling\n")
@time testFullOutput(ffk, Threads.nthreads(), 2.0, verbose)
println("\nParallel: adaptive resampling, essThreshold = 0.75\n")
@time testFullOutput(ffk, Threads.nthreads(), 0.75, verbose)
println("\nParallel: adaptive resampling, essThreshold = 0.25\n")
@time testFullOutput(ffk, Threads.nthreads(), 0.25, verbose)

@time @testset "Finite FK: csmc" begin
  testcsmc(1, 2.0)
  testcsmc(1, 0.5)
  testcsmc(Threads.nthreads(), 2.0)
  testcsmc(Threads.nthreads(), 0.5)
end
