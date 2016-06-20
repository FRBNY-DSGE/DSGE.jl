"""
```
hessizero{T<:AbstractFloat}(fcn::Function, x::Vector{T};
                            check_neg_diag::Bool=false,
                            verbose::Symbol=:none,
                            distr::Bool=true)
```

Compute Hessian of function `fcn` evaluated at `x`.

### Arguments
- `check_neg_diag`: Throw an error if any negative diagonal elements are detected.
- `verbose`: Print verbose output
- `distr`: Use available parallel workers to increase performance.
"""
function hessizero{T<:AbstractFloat}(fcn::Function,
                                    x::Vector{T};
                                    check_neg_diag::Bool=false,
                                    verbose::Symbol=:none,
                                    distr::Bool=true)
    n_para = length(x)
    hessian  = zeros(n_para, n_para)

    # Compute diagonal elements first
    if distr && nworkers() > 1
        diag_elements = @sync @parallel (hcat) for i = 1:n_para
            hess_diag_element(fcn, x, i; check_neg_diag=check_neg_diag, verbose=verbose)
        end
        for i = 1:n_para
            hessian[i, i] = diag_elements[i]
        end
    else
        for i=1:n_para
            hessian[i,i] = hess_diag_element(fcn, x, i; check_neg_diag=check_neg_diag,
                                             verbose=verbose)
        end
    end

    # Now compute off-diagonal elements
    # Make sure that correlations are between -1 and 1
    # invalid_corr indexes elements that are invalid
    invalid_corr = Dict{Tuple{Int,Int}, Float64}()

    # Build indices to iterate over
    n_off_diag_els = Int(n_para*(n_para-1)/2)
    off_diag_inds = Vector{Tuple{Int,Int}}(n_off_diag_els)
    k=1
    for i=1:(n_para-1), j=(i+1):n_para
        off_diag_inds[k] = (i,j)
        k = k+1
    end

    # Iterate over off diag elements
    if distr
        off_diag_out = @sync @parallel (hcat) for (i,j) in off_diag_inds
            σ_xσ_y = sqrt(abs(hessian[i, i]*hessian[j, j]))
            hess_offdiag_element(fcn, x, i, j, σ_xσ_y; verbose=verbose)
        end
        # Ensure off_diag_out is array
        off_diag_out = hcat(off_diag_out)
    else
        off_diag_out = Array{Tuple{T, T},1}(n_off_diag_els)
        for (k,(i,j)) in enumerate(off_diag_inds)
            σ_xσ_y = sqrt(abs(hessian[i, i]*hessian[j, j]))
            off_diag_out[k] = hess_offdiag_element(fcn, x, i, j, σ_xσ_y; verbose=verbose)
        end
    end

    # Fill in values
    for k=1:n_off_diag_els
        (i,j) = off_diag_inds[k]
        (value, ρ_xy) = off_diag_out[k]

        hessian[i,j] = value
        hessian[j,i] = value

        if ρ_xy < -1 || 1 < ρ_xy
            invalid_corr[(i, j)] = ρ_xy
        end
    end

    has_errors = false
    if !isempty(invalid_corr)
        println("Errors: $invalid_corr")
        has_errors = true
    end

    return hessian, has_errors
end

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
