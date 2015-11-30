"""
```
optimize!(model::AbstractDSGEModel, data::Matrix;
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
function optimize!(model::AbstractDSGEModel,
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
        else
            error("Method ",method," is not supported.")
        end
    
        # Inputs to optimization
        H0             = 1e-4 * eye(num_parameters_free(model))
        para_free_inds = find([!θ.fixed for θ in model.parameters])
        x_model        = toreal(model.parameters)
        x_opt          = x_model[para_free_inds]

        function f_opt(x_opt)
            x_model[para_free_inds] = x_opt
            tomodel!(model,x_model)
            return -posterior(model, data; catch_errors=true)[:post]
        end

        rng = model.rng
    
        out, H_ = optimizer(f_opt, x_opt, H0; 
            xtol=xtol, ftol=ftol, grtol=grtol, iterations=iterations,
            store_trace=store_trace, show_trace=show_trace, extended_trace=extended_trace,
            verbose=verbose, rng=rng)

        x_model[para_free_inds] = out.minimum
        tomodel!(model, x_model)

        # Match original dimensions
        out.minimum = x_model

        H = zeros(num_parameters(model), num_parameters(model))
        # Fill in rows/cols of zeros corresponding to location of fixed parameters
        # For each row corresponding to a free parameter, fill in columns corresponding to free
        # parameters. Everything else is 0.
        for (row_free, row_full) in enumerate(para_free_inds)
            H[row_full,para_free_inds] = H_[row_free,:]
        end

        return out, H
end
