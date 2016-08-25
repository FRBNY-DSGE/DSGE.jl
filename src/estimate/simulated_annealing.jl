function simulated_annealing(fcn::Function,
                             x0::Vector,
                             args...;
                             neighbor!::Function  = Optim.default_neighbor!,
                             xtol::Real           = 1e-32, # default from Optim.jl
                             ftol::Float64        = 1e-14, # Default from csminwel
                             iterations::Int      = 1000,
                             store_trace::Bool    = false,
                             show_trace::Bool     = false,
                             extended_trace::Bool = false,
                             kwargs...)

#    println("$(neighbor!)")
    optim_options = OptimizationOptions(x_tol       = xtol,
                                        f_tol       = ftol,
                                        iterations  = iterations,
                                        store_trace = store_trace,
                                        show_trace  = show_trace,
                                        extended_trace = extended_trace)

    optimize(fcn, x0, SimulatedAnnealing(neighbor! = neighbor!), optim_options)
end

# function neighbor_dsge!(x::Array, x_proposal::Array)
#     @assert size(x) == size(x_proposal)

#     for i in 1:length(x)
#         r = randn() * 0.01
#         @inbounds x_proposal[i] = x[i] + r
#     end

    
# end