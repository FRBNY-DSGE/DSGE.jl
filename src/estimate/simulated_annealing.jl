function simulated_annealing(fcn::Function,
                             x0::Array,
                             args...;
                             neighbor!::Function   = Optim.default_neighbor!,
                             xtol::Real            = 1e-32, # default from Optim.jl
                             ftol::Float64         = 1e-14, # Default from csminwel
                             iterations::Int       = 1000,
                             store_trace::Bool     = false,
                             show_trace::Bool      = false,
                             extended_trace::Bool  = false,
                             temperature::Function = Optim.log_temperature,
                             kwargs...)
 
    Optim.optimize(fcn, x0, 
                   method = SimulatedAnnealing(neighbor! = neighbor!,temperature = temperature),
                   iterations = iterations, store_trace = store_trace, show_trace = show_trace,
                   extended_trace = extended_trace), nothing

end


function log_temperature(t::Real; initial_temperature::Real = 1.0)
    return initial_temperature/log(1+t)
end

function exponential_temperature(t::Real; initial_temperature::Real = 1.0, α::Real = .99)
    if !( 0.0 < α < 1.0) 
        error("α must be in (0,1)")
    end
    return initial_temperature*α^t
end

function linear_temperature(t::Real;  initial_temperature::Real = 1, β = .005)
    return max(initial_temperature - t*β, 0.0)
end