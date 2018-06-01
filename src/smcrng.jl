using Compat.Random

# import RandomNumbers.Random123.Threefry4x
#
# const SMCRNG = Threefry4x{UInt64, 20}
#
# let
#   const engines::Vector{SMCRNG} = Vector{SMCRNG}(0)
#   global function getSMCRNG()
#     @inbounds return engines[Threads.threadid()]
#   end
#   global function setSMCRNGs(v::Int64)
#     seed = tuple(0, 0, 0, v)
#     Threads.@threads for i = 1:length(engines)
#       @inbounds engines[Threads.threadid()] =
#         Threefry4x(seed .+ tuple(0, 0, 0, i*1125899906842624))
#     end
#   end
#   ## happens at runtime to avoid false sharing
#   global function initializeSMCRNGs()
#     engines = Vector{SMCRNG}(Threads.nthreads())
#     Threads.@threads for i = 1:length(engines)
#       @inbounds engines[Threads.threadid()] = Threefry4x()
#     end
#   end
# end

const SMCRNG = MersenneTwister

let
  engines::Vector{SMCRNG} = Vector{SMCRNG}(undef, 0)
  global function getSMCRNG()
    @inbounds return engines[Threads.threadid()]
  end
  global function setSMCRNGs(v::Int64)
    Threads.@threads for i = 1:length(engines)
      @inbounds engines[Threads.threadid()] =
        MersenneTwister(v+i*1125899906842624)
    end
  end
  ## happens at runtime to avoid false sharing
  global function initializeSMCRNGs()
    engines = Vector{SMCRNG}(undef, Threads.nthreads())
    Threads.@threads for i = 1:length(engines)
      @inbounds engines[Threads.threadid()] =
        MersenneTwister(Random.make_seed())
    end
    setSMCRNGs(0)
  end
end
