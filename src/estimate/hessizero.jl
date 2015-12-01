function hessizero{T<:AbstractFloat}(fcn::Function, 
                                    x::Vector{T}; 
                                    check_neg_diag::Bool=false,
                                    verbose::Symbol=:none,
                                    distr::Bool=true)
    n_para = length(x)
    hessian  = zeros(n_para, n_para)

    # Compute diagonal elements first
    if distr
        diag_elements = @sync @parallel (hcat) for i = 1:n_para
            hess_diag_element(fcn, x, i; check_neg_diag=check_neg_diag, verbose=verbose)
        end
        for i = 1:n_para
            hessian[i, i] = diag_elements[i]
        end
    else
        for i=1:n_para
            hessian[i,i] = hess_diag_element(fcn, x, i; check_neg_diag=check_neg_diag, verbose=verbose) 
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
