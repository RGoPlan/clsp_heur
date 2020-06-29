module Heuristics

include("Heuristics/relax_and_fix.jl")
include("Heuristics/fix_and_optimize.jl")

export relax_and_fix, fix_and_optimize

end