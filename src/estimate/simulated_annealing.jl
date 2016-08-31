function simulated_annealing(fcn::Function,
                             x0::Array,
                             args...;
                             neighbor!::Function  = Optim.default_neighbor!,
                             xtol::Real           = 1e-32, # default from Optim.jl
                             ftol::Float64        = 1e-14, # Default from csminwel
                             iterations::Int      = 1000,
                             store_trace::Bool    = false,
                             show_trace::Bool     = false,
                             extended_trace::Bool = false,
                             kwargs...)

    Optim.optimize(fcn, x0, method = SimulatedAnnealing(neighbor! = neighbor!),
                   iterations = iterations, store_trace = store_trace, show_trace = show_trace,
                   extended_trace = extended_trace), nothing

end

