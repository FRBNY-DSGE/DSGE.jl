# Compute diag element
function hess_diag_element!{T<:AbstractFloat}(model::AbstractDSGEModel, 
                                               x::Vector{T}, 
                                               YY::Matrix{T}, 
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

        fx  = posterior!(model, x, YY)
        fdx = posterior!(model, paradx, YY)
        fdy = posterior!(model, parady, YY)

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
function hess_offdiag_element!{T<:AbstractFloat}(model::AbstractDSGEModel, 
                                                  x::Vector{T}, 
                                                  YY::Matrix{T}, 
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

        fx    = posterior!(model, x, YY)
        fdx   = posterior!(model, paradx, YY)
        fdy   = posterior!(model, parady, YY)
        fdxdy = posterior!(model, paradxdy, YY)

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

# Compute Hessian of posterior function evaluated at x
function hessian!{T<:AbstractFloat}(model::AbstractDSGEModel, x::Vector{T}, YY::Matrix{T}; verbose::Bool = false)

    update!(model, x)

    ## index of free parameters
    para_free      = [!θ.fixed for θ in model.parameters]
    para_free_inds = find(para_free)
    num_para_free  = length(para_free_inds)

    num_para = length(x)
    hessian  = zeros(num_para, num_para)

    # Compute diagonal elements first
    for row = para_free_inds'
        hessian[row, row] = hess_diag_element!(model, x, YY, row; verbose=verbose)
    end

    # Now compute off-diagonal elements
    # Make sure that correlations are between -1 and 1
    # errorij contains the index of elements that are invalid
    #TODO this can just be a matrix
    errorij = Dict{Tuple{Int64}, Float64}()

    for i = 1:(num_para_free-1)
        row = para_free_inds[i]
        for j = (i+1):num_para_free
            col = para_free_inds[j]

            σ_xσ_y = sqrt(hessian[row, row]*hessian[col, col])
            (value, ρ_xy) = hess_offdiag_element!(model, x, YY, row, col, σ_xσ_y; verbose=verbose)

            hessian[row, col] = value
            hessian[col, row] = value

            # if not null
            if ρ_xy < -1 || ρ_xy > 1
                errorij[(row, col)] = ρ_xy
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
