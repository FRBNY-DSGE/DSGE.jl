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
This routine alternates between LBFGS and simulated annealing.
LBFGS is good for quickly moving in the correct direction, while simulated annealing
is used to escape local minima.
"""
function combined_optimizer(fcn::Function,
                            x0::Vector,
                            args...;
                            xtol::Real            = 1e-32, # default from Optim.jl
                            ftol::Float64         = 1e-14, # default from csminwel
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

    cycle      = 0
    converged  = false
    x_opt      = x0
    use_sa_out = false
    f_opt      = fcn(x_opt)
    out_sa     = nothing
    out_lbfgs  = nothing

    while cycle < max_cycles && converged == false

        # first, run LBFGS
        println(verbose, :low, "Running L-BFGS...")
        out_lbfgs = Optim.optimize(fcn, x_opt, lbfgs(),
                   Optim.Options(autodiff=autodiff, g_tol = grtol, f_tol = ftol, x_tol = xtol,
                   iterations = iterations, store_trace = store_trace, show_trace = show_trace,
                   extended_trace = extended_trace))

        # store relevant information from the optimizer
        minimum_lbfgs   = out_lbfgs.minimum
        minimizer_lbfgs = out_lbfgs.minimizer

        # second, run simulated annealing with the trace disabled
        println(verbose, :low, "Running simulated annealing...")

        # simulated annealing gets more iterations, no trace printout
        sa_iterations = min(iterations * 10, 250)
        out_sa = Optim.optimize(fcn, x_opt,
                                method = SimulatedAnnealing(neighbor!   = neighbor!,
                                                            temperature = temperature),
                                iterations = sa_iterations, store_trace = store_trace,
                                show_trace = false,
                                extended_trace = extended_trace)

        minimum_sa   = out_sa.minimum
        minimizer_sa = out_sa.minimizer

        # reject a move to a worse spot
        use_sa_out = minimum_sa < minimum_lbfgs
        x_opt      = use_sa_out ? minimizer_sa : minimizer_lbfgs
        cycle_best = use_sa_out ? minimum_sa : minimum_lbfgs

        # check convergence
        rel_diff  = abs(f_opt - cycle_best)/abs(f_opt)
        converged = rel_diff < ftol
        cycle    += 1

        println(verbose, :low, "cycle: $cycle")
        println(verbose, :low, "relative function difference: $(round(rel_diff,5))")

        f_opt = cycle_best
    end

    println("optimization complete. cycles: $cycle")
    out = use_sa_out ? out_sa : out_lbfgs
    return out
end
