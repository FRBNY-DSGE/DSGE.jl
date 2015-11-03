# Compute diag element
function hess_diag_element!{T<:AbstractFloat}(fcn::Function,
                                              x::Vector{T}, 
                                              i::Int; 
                                              ndx::Int=6,
                                              verbose::Bool = false)
    # Setup
    num_para = length(x)
    dxscale  = ones(num_para, 1)
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)

    # Computation
    if verbose
        println("Hessian element: ($i, $i)")
    end

    # Diagonal element computation
    for k = 1:ndx
        paradx    = copy(x)
        parady    = copy(x)
        paradx[i] = paradx[i] + dx[k]*dxscale[i]
        parady[i] = parady[i] - dx[k]*dxscale[i]

        fx  = fcn(x)
        fdx = fcn(paradx)
        fdy = fcn(parady)

        hessdiag[k]  = -(2fx - fdx - fdy) / (dx[k]*dxscale[i])^2
    end

    value = -(hessdiag[3]+hessdiag[4])/2

    if value < 0
        error("Negative diagonal in Hessian")
    end
    if verbose
        println("Value used: $value")
    end

    return value
end

# Compute off diag element
function hess_offdiag_element!{T<:AbstractFloat}(fcn::Function,
                                                 x::Vector{T}, 
                                                 i::Int, 
                                                 j::Int,
                                                 σ_xσ_y::T;
                                                 ndx::Int=6,
                                                 verbose::Bool = false)
    # Setup
    num_para = length(x)
    dxscale  = ones(num_para, 1)
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)

    # Computation
    if verbose
        println("Hessian element: ($i, $j)")
    end

    for k = 1:ndx
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

    if verbose
        println("Values: $(-hessdiag)")
    end

    value = -(hessdiag[3]+hessdiag[4])/2

    if value == 0 || σ_xσ_y == 0
        ρ_xy = 0
    else
        ρ_xy = value / σ_xσ_y
    end

    if ρ_xy < -1 || 1 < ρ_xy
        value = 0
    end

    if verbose
        println("Value used: $value")
        println("Correlation: $ρ_xy")
    end

    return value, ρ_xy
end

function hessian!{T<:AbstractFloat}(fcn::Function, 
                                    x::Vector{T}; 
                                    verbose::Bool = false)
    num_para = length(x)
    hessian  = zeros(num_para, num_para)

    # Compute diagonal elements first
    diag_elements = @sync @parallel (hcat) for i = 1:num_para
        hess_diag_element!(fcn, x, i; verbose=verbose)
    end
    for i = 1:num_para
        hessian[i, i] = diag_elements[i]
    end

    # Now compute off-diagonal elements
    # Make sure that correlations are between -1 and 1
    # invalid_corr indexes elements that are invalid
    invalid_corr = Dict{Tuple{Int,Int}, Float64}()

    # Build indices to iterate over
    num_off_diag_els = Int(num_para*(num_para-1)/2)
    off_diag_inds = Vector{Tuple{Int,Int}}(num_off_diag_els)
    k=1
    for i=1:(num_para-1), j=(i+1):num_para
        off_diag_inds[k] = (i,j)
        k = k+1
    end

    # Iterate over off diag elements
    off_diag_out = @sync @parallel (hcat) for (i,j) in off_diag_inds
        σ_xσ_y = sqrt(hessian[i, i]*hessian[j, j])
        hess_offdiag_element!(fcn, x, i, j, σ_xσ_y; verbose=verbose)
    end
    for k=1:num_off_diag_els
        (i,j) = off_diag_inds[k]
        (value, ρ_xy) = off_diag_out[k]

        hessian[i,j] = value
        hessian[j,i] = value

        if ρ_xy < -1 || 1 < ρ_xy
            invalid_corr[(i, j)] = ρ_xy
        end
    end

    has_errors= false
    if !isempty(invalid_corr)
        println("Errors: $invalid_corr")
        has_errors = true
    end

    return hessian, has_errors 
end

# Compute Hessian of posterior function evaluated at x
function hessian!{T<:AbstractFloat}(model::AbstractDSGEModel, 
                                    x::Vector{T}, 
                                    YY::Matrix{T}; 
                                    verbose::Bool = false)
    update!(model, x)

    # Index of free parameters
    para_free      = [!θ.fixed for θ in model.parameters]
    para_free_inds = find(para_free)

    # Compute hessian only for freem parameters with indices less than max. Useful for
    # testing purposes.
    max_free_ind = max_hessian_free_params(model)
    if max_free_ind < maximum(para_free_inds)
        para_free_inds = para_free_inds[1:max_free_ind]
    end

    num_para = length(x)
    hessian  = zeros(num_para, num_para)

    # x_hessian is the vector of free params
    # x_model is the vector of all params
    x_model = copy(x)
    x_hessian = x_model[para_free_inds]
    function f_hessian(x_hessian)
        x_model[para_free_inds] = x_hessian
        return posterior!(model, x_model, YY)
    end

    hessian_free, has_errors = hessian!(f_hessian, x_hessian; verbose=verbose)

    # Fill in rows/cols of zeros corresponding to location of fixed parameters
    # For each row corresponding to a free parameter, fill in columns corresponding to free
    # parameters. Everything else is 0.
    for (row_free, row_full) in enumerate(para_free_inds)
        hessian[row_full,para_free_inds] = hessian_free[row_free,:]
    end

    return hessian, has_errors
end
