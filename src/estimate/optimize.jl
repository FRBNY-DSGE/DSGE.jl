using Debug
# """
# ```
# optimize!(m::AbstractModel, data::Matrix;
#           method::Symbol       = :csminwel,
#           xtol::Real           = 1e-32,  # default from Optim.jl
#           ftol::Float64        = 1e-14,  # Default from csminwel
#           grtol::Real          = 1e-8,   # default from Optim.jl
#           iterations::Int      = 1000,
#           store_trace::Bool    = false,
#           show_trace::Bool     = false,
#           extended_trace::Bool = false,
#           step_size::Float64   = .01,
#           verbose::Symbol      = :none)
# ```

# Wrapper function to send a model to csminwel (or another optimization routine).
# """
@debug function optimize!(m::AbstractModel,
                   data::Matrix;
                   method::Symbol       = :csminwel,
                   xtol::Real           = 1e-32,  # default from Optim.jl
                   ftol::Float64        = 1e-14,  # Default from csminwel
                   grtol::Real          = 1e-8,   # default from Optim.jl
                   iterations::Int      = 1000,
                   store_trace::Bool    = false,
                   show_trace::Bool     = false,
                   extended_trace::Bool = false,
                   step_size::Float64   = .01,
                   verbose::Symbol      = :none)

    # For now, only csminwel should be used
    optimizer = if method == :csminwel
        csminwel
    elseif method == :simulated_annealing
        simulated_annealing
    else
        error("Method ",method," is not supported.")
    end

    # Inputs to optimization
    H0             = 1e-4 * eye(n_parameters_free(m))
    para_free_inds = find([!θ.fixed for θ in m.parameters])
    x_model        = transform_to_real_line(m.parameters)
    x_opt          = x_model[para_free_inds]

    function f_opt(x_opt)
        x_model[para_free_inds] = x_opt
        transform_to_model_space!(m,x_model)
        return -posterior(m, data; catch_errors=true)[:post]
    end


    function neighbor_dsge!(x, x_proposal)
        # This function computes a proposal "next step" during simulated annealing.
        # Inputs:
        # - `x`: current position (of non-fixed states)
        # - `x_proposal`: proposed next position (of non-fixed states). 
        #                 (passed in for pre-allocation purposes)
        # Outputs:
        # - `x_proposal`

        @assert size(x) == size(x_proposal)

        T = eltype(x)
        npara = length(x)

        # Convert x_proposal to model space and expand to full parameter vector
        x_all = T[p.value for p in m.parameters]  # to get fixed values
        x_all[para_free_inds] = x                 # this is from real line

        x_all_model = transform_to_model_space(m.parameters, x_all)
        x_proposal_all = similar(x_all_model)

        success = false
        count = 1
        while !success

            # take a step in model space
            for i in para_free_inds
                prior_var = moments(get(m.parameters[i].prior))[2]
                proposal_in_bounds = false
                proposal = x_all_model[i]
                # draw a new parameter value, and redraw if out of bounds
                while !proposal_in_bounds
                    r = rand([-1 1]) * rand()
                    proposal = x_all_model[i] + (r * step_size * prior_var)
                    if m.parameters[i].valuebounds[1] < proposal &&
                        m.parameters[i].valuebounds[2] > proposal
                        proposal_in_bounds = true
                    end
                end
                @inbounds x_proposal_all[i] = proposal
            
            end

            # check that model can be solved
            try
                update!(m, x_proposal_all)
                #println("trying to solve model in neighbor")
                solve(m)
                #println("done solving model in neighbor")
                x_proposal_all = transform_to_real_line(m.parameters, x_proposal_all)
                success = true
            end
            
        end

        x_proposal[1:end] = x_proposal_all[para_free_inds]

       return
    end

    rng = m.rng
    temperature = get_setting(m, :simulated_annealing_temperature)

    if method == :simulated_annealing
        out, H_ = optimizer(f_opt, x_opt, H0;
                        xtol = xtol, ftol = ftol, grtol = grtol, iterations = iterations, step_size = step_size,
                        store_trace = store_trace, show_trace = show_trace, extended_trace = extended_trace,
                        neighbor! = neighbor_dsge!, verbose = verbose, rng = rng, temperature = temperature)
    elseif method == :csminwel
        out, H_ = optimizer(f_opt, x_opt, H0;
                        xtol = xtol, ftol = ftol, grtol = grtol, iterations = iterations,
                        store_trace = store_trace, show_trace = show_trace, extended_trace = extended_trace,
                        verbose = verbose, rng = rng)
    end

    x_model[para_free_inds] = out.minimizer
    transform_to_model_space!(m, x_model)

    # Match original dimensions
    out.minimizer = x_model

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
