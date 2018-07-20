__precompile__()

module SequentialMonteCarlo

import Compat.undef
import Compat.Statistics.mean

include("smcrng.jl")
include("structures.jl")
include("categorical.jl")
include("common.jl")
include("serial.jl")
include("parallel.jl")
include("interface.jl")

function __init__()
  initializeSMCRNGs()
end

export smc!, csmc!, SMCModel, SMCIO, SMCRNG, getSMCRNG, setSMCRNGs

end # module
