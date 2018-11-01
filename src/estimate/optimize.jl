"""
```
type optimization_result{T}
    minimizer::Vector{T}
    minimum::T
    converged::Bool
    iterations::Int

```
Container type for various optimization outputs
"""
mutable struct optimization_result{T}
    minimizer::Vector{T}
    minimum::T
    converged::Bool
    iterations::Int
end

"""
```
optimize!(m::AbstractModel, data::Matrix;
          method::Symbol       = :csminwel,
          xtol::Real           = 1e-32,  # default from Optim.jl
          ftol::Float64        = 1e-14,  # Default from csminwel
          grtol::Real          = 1e-8,   # default from Optim.jl
          iterations::Int      = 1000,
          store_trace::Bool    = false,
          show_trace::Bool     = false,
          extended_trace::Bool = false,
          verbose::Symbol      = :none)
```

Wrapper function to send a model to csminwel (or another optimization routine).
"""
function optimize!(m::AbstractModel,
                   data::Matrix;
                   method::Symbol       = :csminwel,
                   xtol::Real           = 1e-32,  # default from Optim.jl
                   ftol::Float64        = 1e-14,  # Default from csminwel
                   grtol::Real          = 1e-8,   # default from Optim.jl
                   iterations::Int      = 1000,
                   store_trace::Bool    = false,
                   show_trace::Bool     = false,
                   extended_trace::Bool = false,
                   mle::Bool            = false, # default from estimate.jl
                   step_size::Float64   = .01,
                   verbose::Symbol      = :none)

    ########################################################################################
    ### Step 1: Setup
    ########################################################################################


    # For now, only csminwel should be used
    optimizer = if method == :csminwel
        csminwel
    elseif method == :simulated_annealing
        simulated_annealing
    elseif method == :nelder_mead
        nelder_mead
    elseif method == :combined_optimizer
        combined_optimizer
    elseif method == :lbfgs
        lbfgs
    else
        error("Method ", method, " is not supported.")
    end

    # Inputs to optimization
    H0             = 1e-4 * eye(n_parameters_free(m))
    para_free_inds = find([!θ.fixed for θ in m.parameters])
    x_model        = transform_to_real_line(m.parameters)
    x_opt          = x_model[para_free_inds]

    ########################################################################################
    ### Step 2: Initialize f_opt
    ########################################################################################

    function f_opt(x_opt)
        try
            x_model[para_free_inds] = x_opt
            transform_to_model_space!(m,x_model)
        catch
            return Inf
        end
        if mle
            out = -likelihood(m, data; catch_errors = true)
        else
            out = -posterior(m, data; catch_errors=true)
        end
        out = !isnan(out) ? out : Inf
        return out
    end

    ########################################################################################
    ### Step 3: Optimizer-specific setup, call optimizer
    ########################################################################################

    # variables used across several optimizers
    rng = m.rng
    temperature = get_setting(m, :simulated_annealing_temperature)
    max_cycles = get_setting(m, :combined_optimizer_max_cycles)
    block_frac = get_setting(m, :simulated_annealing_block_proportion)
    H_ = nothing

    function neighbor_dsge!(x, x_proposal)
        # This function computes a proposal "next step" during simulated annealing.
        # Inputs:
        # - `x`: current position (of non-fixed states)
        # - `x_proposal`: proposed next position (of non-fixed states).
        #                 (passed in for pre-allocation purposes)
        # Outputs:
        # - `x_proposal`

        T = eltype(x)
        npara = length(x)
        subset_inds = []
        while length(subset_inds) == 0
            subset_inds = randsubseq(para_free_inds,block_frac)
        end

        # Convert x_proposal to model space and expand to full parameter vector
        x_all = T[p.value for p in m.parameters]  # to get fixed values
        x_all[para_free_inds] = x                 # this is from real line

        x_all_model = transform_to_model_space(m.parameters, x_all)
        x_proposal_all = copy(x_all_model)

        success = false
        while !success
            # take a step in model space
            for i in subset_inds
                prior_var = moments(m.parameters[i])[2]#moments(get(m.parameters[i].prior))[2]
                proposal_in_bounds = false
                proposal = x_all_model[i]
                lower = m.parameters[i].valuebounds[1]
                upper = m.parameters[i].valuebounds[2]
                # draw a new parameter value, and redraw if out of bounds
                while !proposal_in_bounds
                    r = rand([-1 1]) * rand()
                    proposal = x_all_model[i] + (r * step_size * prior_var)
                    if lower < proposal < upper
                        proposal_in_bounds = true
                    end
                end
                @inbounds x_proposal_all[i] = proposal
            end

            # check that model can be solved
            try
                update!(m, x_proposal_all)
                solve(m)
                x_proposal_all = transform_to_real_line(m.parameters, x_proposal_all)
                success = true
            catch ex
                if !(typeof(ex) in [DomainError, ParamBoundsError, GensysError])
                    rethrow(ex)
                end
            end
        end

        x_proposal[1:end] = x_proposal_all[para_free_inds]

        return
    end

    rng = m.rng
    temperature = get_setting(m, :simulated_annealing_temperature)
    if method == :simulated_annealing
       opt_result = optimizer(f_opt, x_opt;
                        iterations = iterations, step_size = step_size,
                        store_trace = store_trace, show_trace = show_trace, extended_trace = extended_trace,
                        neighbor! = neighbor_dsge!, verbose = verbose, rng = rng, temperature = temperature)
       converged = opt_result.iteration_converged
       out = optimization_result(opt_result.minimizer, opt_result.minimum, converged, opt_result.iterations)

    elseif method == :nelder_mead
        opt_result = optimizer(f_opt, x_opt;
                               iterations = iterations,
                               store_trace = store_trace, show_trace = show_trace,
                               extended_trace = extended_trace, verbose = verbose, rng = rng)
       converged = opt_result.iteration_converged
       out = optimization_result(opt_result.minimizer, opt_result.minimum, converged, opt_result.iterations)

    elseif method == :csminwel
        opt_result, H_ = optimizer(f_opt, x_opt, H0;
                        xtol = xtol, ftol = ftol, grtol = grtol, iterations = iterations,
                        store_trace = store_trace, show_trace = show_trace, extended_trace = extended_trace,
                        verbose = verbose, rng = rng)
        converged = opt_result.g_converged || opt_result.f_converged #|| opt_result.x_converged
        out = optimization_result(opt_result.minimizer, opt_result.minimum, converged, opt_result.iterations)

    elseif method == :lbfgs
        opt_result = optimizer(f_opt, x_opt;
                        xtol = xtol, ftol = ftol, grtol = grtol, iterations = iterations,
                        store_trace = store_trace, show_trace = show_trace, extended_trace = extended_trace,
                        verbose = verbose, rng = rng)
        converged = opt_result.g_converged || opt_result.f_converged #|| opt_result.x_converged
        out = optimization_result(opt_result.minimizer, opt_result.minimum, converged, opt_result.iterations)

    elseif method == :combined_optimizer
        opt_result = optimizer(f_opt, x_opt;
                        xtol = xtol, ftol = ftol, grtol = grtol, iterations = iterations, step_size = step_size,
                        store_trace = store_trace, show_trace = show_trace, extended_trace = extended_trace,
                        neighbor! = neighbor_dsge!, verbose = verbose, rng = rng, temperature = temperature,
                               max_cycles = max_cycles)
        converged = opt_result.g_converged || opt_result.f_converged || opt_result.x_converged
        converged = opt_result.method == "Simulated Annealing" ? opt_result.iteration_converged : converged
        out = optimization_result(opt_result.minimizer, opt_result.minimum, converged, opt_result.iterations)
    end

    ########################################################################################
    ### Step 4: transform output, populate Hessian
    ########################################################################################

    x_model[para_free_inds] = out.minimizer
    transform_to_model_space!(m, x_model)

    # Match original dimensions
    out.minimizer = map(θ -> θ.value, m.parameters)

    H = zeros(n_parameters(m), n_parameters(m))
    if H_ != nothing

        # Fill in rows/cols of zeros corresponding to location of fixed parameters
        # For each row corresponding to a free parameter, fill in columns corresponding to
        # free parameters. Everything else is 0.
        for (row_free, row_full) in enumerate(para_free_inds)
            H[row_full,para_free_inds] = H_[row_free,:]
        end

    end

    return out, H
end
