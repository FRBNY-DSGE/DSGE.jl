### HOW TO RUN THIS CODE
#
# From any directory:
#
#   using DSGE
#   using .M990
#   model = Model()

module DSGE

export
    # modules
    M990,

    # functions
    solve

include("init/DistributionsExt.jl")
include("init/FinancialFrictionsFunctions.jl")
include("AbstractModel.jl")

include("solve/Gensys.jl")
include("solve/solve.jl")

include("estimate/Kalman.jl")

include("models/M990.jl")

end
