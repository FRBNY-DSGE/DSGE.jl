### HOW TO RUN THIS CODE
#
# From any directory:
#
#   using DSGE
#   using .M990
#   model = Model()

module DSGE

export M990

include("init/DistributionsExt.jl")
include("init/FinancialFrictionsFunctions.jl")
include("AbstractModel.jl")

include("solve/Gensys.jl")
include("solve/Kalman.jl")

include("models/M990.jl")

end
