"""
```
hessian!{T<:AbstractFloat}(m::AbstractModel, x::Vector{T}, data::Matrix{T};
                           verbose::Symbol = :none)
```

Compute Hessian of DSGE posterior function evaluated at x.
"""
function hessian!{T<:AbstractFloat}(m::AbstractModel,
                                    x::Vector{T},
                                    data::Matrix{T};
                                    verbose::Symbol = :none)
    update!(m, x)

    # Index of free parameters
    para_free      = [!θ.fixed for θ in m.parameters]
    para_free_inds = find(para_free)

    # Compute hessian only for freem parameters with indices less than max. Useful for
    # testing purposes.
    max_free_ind = n_hessian_test_params(m)
    if max_free_ind < maximum(para_free_inds)
        para_free_inds = para_free_inds[1:max_free_ind]
    end

    n_para = length(x)
    hessian  = zeros(n_para, n_para)

    # x_hessian is the vector of free params
    # x_model is the vector of all params
    x_model = copy(x)
    x_hessian = x_model[para_free_inds]
    function f_hessian(x_hessian)
        x_model[para_free_inds] = x_hessian
        return -posterior!(m, x_model, data)[:post]
    end

    distr=use_parallel_workers(m)
    hessian_free, has_errors = hessizero(f_hessian, x_hessian;
        check_neg_diag=true, verbose=verbose, distr=distr)

    # Fill in rows/cols of zeros corresponding to location of fixed parameters
    # For each row corresponding to a free parameter, fill in columns corresponding to free
    # parameters. Everything else is 0.
    for (row_free, row_full) in enumerate(para_free_inds)
        hessian[row_full,para_free_inds] = hessian_free[row_free,:]
    end

    return hessian, has_errors
end
