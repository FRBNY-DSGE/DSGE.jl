# Compute diag element
function hess_diag_element{T<:AbstractFloat}(fcn::Function,
                                              x::Vector{T}, 
                                              i::Int; 
                                              ndx::Int=6,
                                              check_neg_diag::Bool=false,
                                              verbose::Symbol=:none)
    # Setup
    n_para = length(x)
    dxscale  = ones(n_para, 1)
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)

    # Computation
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Hessian element: ($i, $i)")
    end

    # Diagonal element computation
    for k = 3:4
        paradx    = copy(x)
        parady    = copy(x)
        paradx[i] = paradx[i] + dx[k]*dxscale[i]
        parady[i] = parady[i] - dx[k]*dxscale[i]

        fx  = fcn(x)
        fdx = fcn(paradx)
        fdy = fcn(parady)

        hessdiag[k]  = -(2fx - fdx - fdy) / (dx[k]*dxscale[i])^2
    end

    if VERBOSITY[verbose] >= VERBOSITY[:high]
        println("Values: $(hessdiag)")
    end

    value = (hessdiag[3]+hessdiag[4])/2

    if check_neg_diag && value < 0
        error("Negative diagonal in Hessian")
    end

    if VERBOSITY[verbose] >= VERBOSITY[:high]
        println("Value used: $value")
    end

    return value
end

# Compute off diag element
function hess_offdiag_element{T<:AbstractFloat}(fcn::Function,
                                                 x::Vector{T}, 
                                                 i::Int, 
                                                 j::Int,
                                                 σ_xσ_y::T;
                                                 ndx::Int=6,
                                                 verbose::Symbol=:none)
    # Setup
    n_para = length(x)
    dxscale  = ones(n_para, 1)
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)

    # Computation
    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Hessian element: ($i, $j)")
    end

    for k = 3:4
        paradx      = copy(x)
        parady      = copy(x)
        paradx[i]   = paradx[i] + dx[k]*dxscale[i]
        parady[j]   = parady[j] - dx[k]*dxscale[j]
        paradxdy    = copy(paradx)
        paradxdy[j] = paradxdy[j] - dx[k]*dxscale[j]

        fx    = fcn(x)
        fdx   = fcn(paradx)
        fdy   = fcn(parady)
        fdxdy = fcn(paradxdy)

        hessdiag[k]  = -(fx - fdx - fdy + fdxdy) / (dx[k]*dx[k]*dxscale[i]*dxscale[j])
    end

    if VERBOSITY[verbose] >= VERBOSITY[:high]
        println("Values: $(hessdiag)")
    end

    value = (hessdiag[3]+hessdiag[4])/2

    if value == 0 || σ_xσ_y == 0
        ρ_xy = 0
    else
        ρ_xy = value / σ_xσ_y
    end

    if ρ_xy < -1 || 1 < ρ_xy
        value = 0
    end

    if VERBOSITY[verbose] >= VERBOSITY[:high]
        println("Value used: $value")
        println("Correlation: $ρ_xy")
    end

    return value, ρ_xy
end

# Compute Hessian of posterior function evaluated at x
function hessian!{T<:AbstractFloat}(m::AbstractModel, 
                                    x::Vector{T}, 
                                    YY::Matrix{T}; 
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
        return -posterior!(m, x_model, YY)[:post]
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
