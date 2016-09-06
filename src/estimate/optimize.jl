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
                   verbose::Symbol      = :none)

    # For now, only csminwel should be used
    optimizer = if method == :csminwel
        csminwel
    elseif method == :simulatedannealing
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
        println("Calling posterior...")
        res= -posterior(m, data; catch_errors=true)[:post]
        println("Done calling posterior")
        return res
    end

    
    function neighbor_dsge!(x, x_proposal; cc = 0.01)
        # This function computes a proposal "next step" during simulated annealing.
        # Inputs:
        # - `x`: current position (of non-fixed states)
        # - `x_proposal`: proposed next position (of non-fixed states).
        # Outputs:
        # - `x_proposal`
        
        @assert size(x) == size(x_proposal)

        T = eltype(x)
        npara = length(x)

        # get covariance matrix from prior
        prior_draws = rand_prior(m)
        prior_cov   = cov(prior_draws)
        
        # Convert x_proposal to model space and expand to full parameter vector
        x_all = T[p.value for p in m.parameters]  # to get fixed values
        x_all[para_free_inds] = x                 # this is from real line
            
        x_all_model = transform_to_model_space(m.parameters, x_all)
        x_proposal_all = similar(x_all_model)
        
        success = false
        while !success

            # take a step in model space
            for i in para_free_inds
                r = rand([-1 1]) * rand()
                @inbounds x_proposal_all[i] = x_all_model[i] + (r * cc * prior_cov[i,i])
            end
                   
            # check that parameters are inbounds, model can be solved,
            # and parameters can be transformed to the real line.
            #try
                update!(m, x_proposal_all)
                println("trying to solve model in neighbor")
                solve(m)
                println("done solving model in neighbor")
                x_proposal_all = transform_to_real_line(m.parameters, x_proposal_all)
                success = true
            #catch err
                #success = false
                # rethrow(err)
                #warn("There was a $(typeof(err))")
            #end
        end
    
        # extract free inds
        # fixed_count = 0
        # for i in 1:length(para_free_inds)
        #     if !in(i, para_free_inds)
        #         fixed_count += 1
        #     end
        #     x_proposal[i] = x_proposal_all[i+fixed_count]
        # end
        x_proposal[1:end] = x_proposal_all[para_free_inds]    

       return
    end

    rng = m.rng

    out, H_ = optimizer(f_opt, x_opt, H0;
                        xtol = xtol, ftol = ftol, grtol = grtol, iterations = iterations,
                        store_trace = store_trace, show_trace = show_trace, extended_trace = extended_trace,
                        neighbor! = neighbor_dsge!, verbose = verbose, rng = rng)
    
    x_model[para_free_inds] = out.minimum
    transform_to_model_space!(m, x_model)
    
    # Match original dimensions
    out.minimum = x_model

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
