function lbfgs(fcn::Function,
                       x0::Vector,
                       args...;
                       xtol::Real           = 1e-32, # default from Optim.jl
                       ftol::Float64        = 1e-14, # Default from csminwel
                       grtol::Real          = 1e-8,  # default from Optim.jl
                       iterations::Int      = 1000,
                       store_trace::Bool    = false,
                       show_trace::Bool     = false,
                       extended_trace::Bool = false,
                       verbose::Symbol      = :none,
                       rng::AbstractRNG     = MersenneTwister(),
                       autodiff::Bool       = false,
                       kwargs...)

    Optim.optimize(fcn, x0,
                   LBFGS(),
                   Optim.Options(autodiff=autodiff, g_tol = grtol, f_tol = ftol, x_tol = xtol,
                   iterations = iterations, store_trace = store_trace, show_trace = show_trace,
                   extended_trace = extended_trace))
end
