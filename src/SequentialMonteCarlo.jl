VERSION < v"0.7-" && __precompile__()

module SequentialMonteCarlo

import Compat.undef
import Compat.Statistics.mean
using RNGPool

include("structures.jl")
include("categorical.jl")
include("common.jl")
include("serial.jl")
include("parallel.jl")
include("interface.jl")

export smc!, csmc!, SMCModel, SMCIO

end # module
