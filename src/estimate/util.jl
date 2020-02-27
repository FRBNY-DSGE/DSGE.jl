"""
```
function thin_mh_draws(m::AbstractDSGEModel, params::Matrix{Float64}; jstep::Int64 = 1)
```
Allows for a custom thinning factor (jstep) to be specified.
If not, it pulls the jstep from the model
"""
function thin_mh_draws(m::AbstractDSGEModel, params::Matrix{Float64}; jstep::Int64 = 1)
    jstep = jstep == 1 ? m.settings[:forecast_jstep].value : jstep
    n_total_draws, n_params = size(params)

    # Thin as usual if n_total_draws % jstep == 0
    # If it does not evenly divide, then start from the remainder+1-th index
    # and then take a thinned subset from there
    n_new_draws, offset = divrem(n_total_draws, jstep)
    params_thinned = Matrix{Float64}(undef, n_new_draws, n_params)
    params_offset = params[offset+1:end, :]

    for (i, j) in enumerate(1:jstep:n_total_draws)
        params_thinned[i, :] = params_offset[j, :]
    end
    return params_thinned
end

"""
```
function calculate_mode(m::Union{AbstractDSGEModel, AbstractVARModel}, data::Matrix{S},
                        starting_mode::Vector{S}, method::Symbol = :csminwel; mle::Bool = false,
                        verbose::Symbol = :none) where {S <: Real}
```
uses a numerical optimizer to calculate the mode of the posterior distribution
starting from pre-specified parameters. This function automatically saves
the mode found to an HDF5 file.
"""
function calculate_mode(m::Union{AbstractDSGEModel, AbstractVARModel}, data::Matrix{S},
                        starting_params::Vector{S}, method::Symbol = :csminwel;
                        mle::Bool = false, get_all_results::Bool = false,
                        save_results::Bool = true, check_neg_diag::Bool = true,
                        verbose::Symbol = :none) where {S <: Real}
    n_iterations = get_setting(m, :optimization_iterations)
    ftol         = get_setting(m, :optimization_ftol)
    xtol         = get_setting(m, :optimization_xtol)
    gtol         = get_setting(m, :optimization_gtol)
    step_size    = get_setting(m, :optimization_step_size)
    max_attempts = get_setting(m, :optimization_attempts)

    DSGE.update!(m, starting_params)

    modal_params  = optimize_step!(m, data, method; mle = mle, n_iterations = n_iterations,
                                   ftol = ftol, xtol = xtol, gtol = gtol, step_size = step_size,
                                   max_attempts = max_attempts, get_all_results = get_all_results,
                                   save_results = save_results, verbose = verbose)
    modal_hessian = hessian_step!(m, data; check_neg_diag = check_neg_diag, verbose = verbose)

    if get_all_results
        return modal_params..., modal_hessian
    else
        return modal_params, modal_hessian
    end
end

"""
```
function optimize_step!(m::Union{AbstractDSGEModel, AbstractVARModel}, data::Matrix{S},
                       method::Symbol = :csminwel; mle::Bool = false, n_iterations::Int = 4,
                       ftol::S = 1e-10, xtol::S = 1e-10, gtol::S = 1e-10, step_size::S = 0.01,
                       max_attempts::Int = 4, verbose::Symbol = :none) where {S <: Real}
```
maximizes the posterior or likelihood of the model object passed using
a numerical optimizer, which by default is csminwel. This function
automatically saves the mode found.
"""
function optimize_step!(m::Union{AbstractDSGEModel, AbstractVARModel}, data::Matrix{S},
                        method::Symbol = :csminwel; mle::Bool = false, n_iterations::Int = 4,
                        ftol::S = 1e-10, xtol::S = 1e-10, gtol::S = 1e-10, step_size::S = 0.01,
                        max_attempts::Int = 4, get_all_results::Bool = false,
                        save_results::Bool = true,
                        verbose::Symbol = :none) where {S <: Real}
    converged = false

    # If the algorithm stops only because we have exceeded the maximum number of
    # iterations, continue improving guess of modal parameters
    total_iterations = 0
    optimization_time = 0
    attempts = 1

    params = map(θ -> θ.value, get_parameters(m)) # define in global scope first
    out = nothing
    H   = nothing

    while !converged
        begin_time = time_ns()
        out, H = optimize!(m, data; method = method,
                           ftol = ftol, grtol = gtol, xtol = xtol,
                           iterations = n_iterations, show_trace =
                           (verbose == :none ? false : true), step_size = step_size,
                           verbose = verbose,
                           mle = mle)

        attempts += 1
        total_iterations += out.iterations
        converged = out.converged || attempts > max_attempts

        end_time = (time_ns() - begin_time)/1e9
        println(verbose, :low, @sprintf "Total iterations completed: %d\n" total_iterations)
        println(verbose, :low,
                @sprintf "Optimization time elapsed: %5.2f\n" optimization_time += end_time)

        # Write params to file after every `n_iterations` iterations
        params = map(θ -> θ.value, get_parameters(m))
        if save_results
            h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
                file["params"] = params
            end
        end
    end

    # write parameters to file one last time so we have the final mode
    if save_results
        h5open(rawpath(m, "estimate", lowercase(string(get_setting(m, :sampling_method))) *
                       "_" * "paramsmode.h5"),"w") do file
                           file["params"] = params
                       end
    end

    if get_all_results
        return params, out, H
    else
        return params
    end
end

"""
```
function hessian_step!(m::Union{AbstractDSGEModel, AbstractVARModel},
                       data::Matrix{S}; verbose::Symbol = :none,
                       check_neg_diag::Bool = true) where {S <: Real}
```
calculates the Hessian at the current calibration of the model object
and checks that it is positive definite.
"""
function hessian_step!(m::Union{AbstractDSGEModel, AbstractVARModel},
                       data::Matrix{S}; verbose::Symbol = :none,
                       check_neg_diag::Bool = true) where {S <: Real}

    params = map(θ -> θ.value, get_parameters(m))
    hessian, _ = hessian!(m, params, data; verbose = verbose, check_neg_diag = check_neg_diag)


    return hessian
end
