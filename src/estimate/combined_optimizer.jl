"""
```
combined_optimizer(fcn::Function,
                       x0::Vector,
                       args...;
                       xtol::Real            = 1e-32, # default from Optim.jl
                       ftol::Float64         = 1e-14, # Default from csminwel
                       grtol::Real           = 1e-8,  # default from Optim.jl
                       iterations::Int       = 1000,
                       store_trace::Bool     = false,
                       show_trace::Bool      = false,
                       extended_trace::Bool  = false,
                       verbose::Symbol       = :none,
                       rng::AbstractRNG      = MersenneTwister(),
                       autodiff::Bool        = false,
                       neighbor!::Function   = Optim.default_neighbor!,
                       temperature::Function = Optim.log_temperature,
                       max_cycles::Int       = 4,
                       kwargs...)

```
This routine alternates between LBFGS and simulated annealing. LBFGS is good for quickly moving in the correct
direction, while simulated annealing is used to escape local minimima.
"""
function combined_optimizer(fcn::Function,
                       x0::Vector,
                       args...;
                       xtol::Real            = 1e-32, # default from Optim.jl
                       ftol::Float64         = 1e-14, # Default from csminwel
                       grtol::Real           = 1e-8,  # default from Optim.jl
                       iterations::Int       = 1000,
                       store_trace::Bool     = false,
                       show_trace::Bool      = false,
                       extended_trace::Bool  = false,
                       verbose::Symbol       = :none,
                       rng::AbstractRNG      = MersenneTwister(),
                       autodiff::Bool        = false,
                       neighbor!::Function   = Optim.default_neighbor!,
                       temperature::Function = Optim.log_temperature,
                       max_cycles::Int       = 4,
                       kwargs...)


    cycle = 0
    converged = false
    x_opt = x0
    use_SA_out = false
    f_opt = fcn(x_opt)
    out_SA = nothing
    out_LBFGS = nothing

    while cycle < max_cycles && converged == false

        # first, run LBFGS
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Running L-BFGS...")
        end
        out_LBFGS = Optim.optimize(fcn, x_opt, LBFGS(),
                   Optim.Options(autodiff=autodiff, g_tol = grtol, f_tol = ftol, x_tol = xtol,
                   iterations = iterations, store_trace = store_trace, show_trace = show_trace,
                   extended_trace = extended_trace))

        # store relevant information from the optimizer
        minimum_LBFGS = out_LBFGS.minimum
        minimizer_LBFGS = out_LBFGS.minimizer

        # second, run simulated annealing with the trace disabled
        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("Running simulated annealing...")
        end

        # simulated annealing gets more iterations, no trace printout
        SA_iterations = min(iterations * 10, 250)
        out_SA = Optim.optimize(fcn, x_opt,
                             method = SimulatedAnnealing(neighbor! = neighbor!,temperature = temperature),
                             iterations = SA_iterations, store_trace = store_trace, show_trace = false,
                             extended_trace = extended_trace)

        minimum_SA = out_SA.minimum
        minimizer_SA = out_SA.minimizer

        # reject a move to a worse spot
        use_SA_out = minimum_SA < minimum_LBFGS
        x_opt =  use_SA_out ? minimizer_SA : minimizer_LBFGS
        cycle_best = use_SA_out ? minimum_SA : minimum_LBFGS

        # check convergence
        rel_diff = abs(f_opt - cycle_best)/abs(f_opt)
        converged = rel_diff < ftol
        cycle += 1

        if VERBOSITY[verbose] >= VERBOSITY[:low]
            println("cycle: $cycle")
            println("relative function difference: $(round(rel_diff,5))")
        end

        f_opt = cycle_best

    end

    println("optimization complete. cycles: $cycle")
    out = use_SA_out ? out_SA : out_LBFGS
    return out

end
