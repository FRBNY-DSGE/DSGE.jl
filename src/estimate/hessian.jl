# Compute diag element
function hessiandiagelement!{T<:AbstractFloat}(model::AbstractDSGEModel, 
                                               x::Vector{T}, 
                                               YY::Matrix{T}, 
                                               i::Int; 
                                               verbose::Bool = false)
    # Setup
    npara    = length(x)
    ndx      = 6
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)
    dxscale  = ones(npara, 1)

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
function hessianoffdiagelement!{T<:AbstractFloat}(model::AbstractDSGEModel, 
                                                  x::Vector{T}, 
                                                  YY::Matrix{T}, 
                                                  i::Int, 
                                                  j::Int,
                                                  σ_xσ_y::T;
                                                  verbose::Bool = false)
    # Setup
    npara    = length(x)
    ndx      = 6
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessdiag = zeros(ndx, 1)
    dxscale  = ones(npara, 1)

    # Computation
    if verbose
        println("Hessian element: ($seli, $selj)")
    end

    for k = 1:ndx
        paradx         = copy(x)
        parady         = copy(x)
        paradx[seli]   = paradx[seli] + dx[k]*dxscale[seli]
        parady[selj]   = parady[selj] - dx[k]*dxscale[selj]
        paradxdy       = copy(paradx)
        paradxdy[selj] = paradxdy[selj] - dx[k]*dxscale[selj]

        fx    = posterior!(model, x, YY)
        fdx   = posterior!(model, paradx, YY)
        fdy   = posterior!(model, parady, YY)
        fdxdy = posterior!(model, paradxdy, YY)

        hessdiag[k]  = -(fx - fdx - fdy + fdxdy) / (dx[k]*dx[k]*dxscale[seli]*dxscale[selj])
    end

    if verbose
        println("Values: $(-hessdiag)")
    end

    val = -(hessdiag[3]+hessdiag[4])/2

    if val == 0 || σ_xσ_y == 0
        ρ_xy = 0
    else
        ρ_xy = val / σ_xσ_y
    end

    if ρ_xy < -1 || ρ_xy > 1
        val = 0
    end

    return value, ρ_xy
end

# Compute Hessian of posterior function evaluated at x (vector)
# if verbose, display error messages, results, etc.
# 11/12/01 translated by Marco DelNegro in matlab from Frank Schorfheide's program in gauss
function hessizero!{T<:AbstractFloat}(model::AbstractDSGEModel, x::Vector{T}, YY::Matrix{T}; verbose::Bool = false)

    update!(model, x)

    ## index of free parameters
    para_free  = [!θ.fixed for θ in model.parameters]
    fpara_free = find(para_free)
    nfree      = length(fpara_free)

    npara    = length(x)
    ndx      = 6
    dx       = exp(-(6:2:(6+(ndx-1)*2))')
    hessian  = zeros(npara, npara)
    hessdiag = zeros(ndx, 1)
    dxscale  = ones(npara, 1)

    # Compute diagonal elements first
    for seli = fpara_free'
        hessian[seli, seli] = hessiandiagelement!(model, x, YY, seli; verbose=verbose)
    end

    # Now compute off-diagonal elements
    # Make sure that correlations are between -1 and 1
    # errorij contains the index of elements that are invalid
    #TODO this can just be a matrix
    errorij = Dict{Tuple{Int64}, Float64}()

    for i = 1:(nfree-1)
        seli = fpara_free[i]
        for j = (i+1):nfree
            selj = fpara_free[j]

            σ_xσ_y = sqrt(hessian[seli, seli]*hessian[selj, selj])
            (val, ρ_xy) = hessianoffdiagelement!(model, x, YY, seli, selj, σ_xσ_y; verbose=verbose)

            hessian[seli, selj] = val
            hessian[selj, seli] = val

            # if not null
            if ρ_xy < -1 || ρ_xy > 1
            errorij[(seli, selj)] = ρ_xy
            end

            if verbose
                println("Value used: $val")
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
