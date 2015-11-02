# Compute diag element
function hess_diag_element!{T<:AbstractFloat}(fcn::Function,
                                              x::Vector{T}, 
                                              i::Int; 
                                              verbose::Bool = false)
    # Setup
    num_para = length(x)
    ndx      = 6
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)
    dxscale  = ones(num_para, 1)

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
        println("Value used: $value\n")
    end

    return value
end

# Compute off diag element
function hess_offdiag_element!{T<:AbstractFloat}(fcn::Function,
                                                 x::Vector{T}, 
                                                 i::Int, 
                                                 j::Int,
                                                 σ_xσ_y::T;
                                                 verbose::Bool = false)
    # Setup
    num_para = length(x)
    ndx      = 6
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)
    dxscale  = ones(num_para, 1)

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

    if ρ_xy < -1 || ρ_xy > 1
        value = 0
    end

    return value, ρ_xy
end
function hessian!{T<:AbstractFloat}(fcn::Function, 
                                    x::Vector{T}; 
                                    verbose::Bool = false)
    num_para = length(x)
    hessian  = zeros(num_para, num_para)

    # Compute diagonal elements first
    for i = 1:num_para
        hessian[i, i] = hess_diag_element!(fcn, x, i; verbose=verbose)
    end

    # Now compute off-diagonal elements
    # Make sure that correlations are between -1 and 1
    # errorij contains the index of elements that are invalid
    #TODO this can just be a matrix
    errorij = Dict{Tuple{Int64}, Float64}()

    for i = 1:(num_para-1)
        for j = (i+1):num_para

            σ_xσ_y = sqrt(hessian[i, i]*hessian[j, j])
            (value, ρ_xy) = hess_offdiag_element!(fcn, x, i, j, σ_xσ_y; verbose=verbose)

            hessian[i, j] = value
            hessian[j, i] = value

            # if not null
            if ρ_xy < -1 || ρ_xy > 1
                errorij[(i, j)] = ρ_xy
            end

            if verbose
                println("Value used: $value")
                println("Correlation: $ρ_xy")
                println("Number of errors: $(length(errorij))\n")
            end
        end
    end

    has_errors= false
    if !isempty(errorij)
        println("Errors: $errorij")
        has_errors = true
    end

    return hessian, has_errors 
end

# Compute Hessian of posterior function evaluated at x
function hessian!{T<:AbstractFloat}(model::AbstractDSGEModel, x::Vector{T}, YY::Matrix{T}; verbose::Bool = false)
    update!(model, x)

    ## index of free parameters
    para_free      = [!θ.fixed for θ in model.parameters]
    para_free_inds = find(para_free)
    num_para_free  = length(para_free_inds)

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
    for i=1:length(para_free_inds)
        row_full = para_free_inds[i]
        row_free = i
        hessian[row_full,para_free_inds] = hessian_free[row_free,:]
    end

    return hessian, has_errors
end
