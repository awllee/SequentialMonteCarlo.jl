using Compat.Test
import Compat.undef

function testapproxequal(a::Vector{Float64}, b::Vector{Float64}, tol::Float64,
  verbose::Bool)
  v::Float64 = maximum(abs.(a-b))
  verbose && (println(a) ; println(b))
  verbose && println("$v < $tol ?")
  @test v < tol
end

function testapproxequal(a::Float64, b::Float64, tol::Float64, verbose::Bool)
  v::Float64 = abs(a-b)
  verbose && println("$v < $tol ?")
  @test v < tol
end

include("finiteFK_test.jl")

@testset "LGModel tests" begin
  @time include("lgModel_test.jl")
end

@testset "MVLGModel tests" begin
  @time include("mvlgModel_test.jl")
end

@testset "SMC Sampler tests" begin
  @time include("smcSampler_test.jl")
end
