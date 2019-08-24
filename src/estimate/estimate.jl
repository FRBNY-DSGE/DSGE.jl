"""
```
estimate(m, data; verbose=:low, proposal_covariance=Matrix(), method=:SMC)
```

Estimate the DSGE parameter posterior distribution.

### Arguments:
- `m::AbstractDSGEModel`: model object

### Optional Arguments:
- `data`: well-formed data as `Matrix` or `DataFrame`. If this is not provided, the `load_data` routine will be executed.

### Keyword Arguments:
- `verbose::Symbol`: The desired frequency of function progress messages printed to standard out.
   - `:none`: No status updates will be reported.
   - `:low`: Status updates will be provided in csminwel and at each block in
     Metropolis-Hastings.
   - `:high`: Status updates provided at each iteration in Metropolis-Hastings.
- `proposal_covariance::Matrix`: Used to test the metropolis_hastings algorithm with a precomputed
  covariance matrix for the proposal distribution. When the Hessian is singular,
  eigenvectors corresponding to zero eigenvectors are not well defined, so eigenvalue
  decomposition can cause problems. Passing a precomputed matrix allows us to ensure that
  the rest of the routine has not broken.
- `method::Symbol`: The method to use when sampling from the posterior distribution. Can
    be either `:MH` for standard Metropolis Hastings Markov Chain Monte Carlo, or `:SMC`
    for Sequential Monte Carlo.
- `mle`: Set to true if parameters should be estimated by maximum likelihood directly.
    If this is set to true, this function will return after estimating parameters.
- `sampling`: Set to false to disable sampling from the posterior.
"""
function estimate(m::AbstractDSGEModel, df::DataFrame;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0, 0),
                  mle::Bool = false,
                  sampling::Bool = true,
                  filestring_addl::Vector{String} = Vector{String}(),
                  continue_intermediate::Bool = false,
                  intermediate_stage_start::Int = 0,
                  intermediate_stage_increment::Int = 10,
                  save_intermediate::Bool = false)
    data = df_to_matrix(m, df)
    estimate(m, data; verbose = verbose, proposal_covariance = proposal_covariance,
             mle = mle, sampling = sampling,
             intermediate_stage_increment = intermediate_stage_increment,
             save_intermediate = save_intermediate)
end

function estimate(m::AbstractDSGEModel;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0, 0),
                  mle::Bool = false,
                  sampling::Bool = true,
                  filestring_addl::Vector{String} = Vector{String}(),
                  continue_intermediate::Bool = false,
                  intermediate_stage_start::Int = 0,
                  intermediate_stage_increment::Int = 10,
		          save_intermediate::Bool = false)
    # Load data
    df = load_data(m; verbose = verbose)
    estimate(m, df; verbose = verbose, proposal_covariance = proposal_covariance,
             mle = mle, sampling = sampling,
             intermediate_stage_increment = intermediate_stage_increment,
	         save_intermediate = save_intermediate)
end

function estimate(m::AbstractDSGEModel, data::AbstractArray;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0,0),
                  mle::Bool = false,
                  sampling::Bool = true,
                  filestring_addl::Vector{String} = Vector{String}(),
                  continue_intermediate::Bool = false,
                  intermediate_stage_start::Int = 0,
                  intermediate_stage_increment::Int = 10,
		          save_intermediate::Bool = false)

    if !(get_setting(m, :sampling_method) in [:SMC, :MH])
        error("method must be :SMC or :MH")
    else
        method = get_setting(m, :sampling_method)
    end

    ########################################################################################
    ### Step 1: Find posterior/likelihood mode (if reoptimizing, run optimization routine)
    ########################################################################################

    # Specify starting mode

    vint = get_setting(m, :data_vintage)
    if reoptimize(m) && method == :MH
        println("Reoptimizing...")

        # Inputs to optimization algorithm
        n_iterations       = get_setting(m, :optimization_iterations)
        ftol               = get_setting(m, :optimization_ftol)
        xtol               = get_setting(m, :optimization_xtol)
        gtol               = get_setting(m, :optimization_gtol)
        step_size          = get_setting(m, :optimization_step_size)
        converged          = false

        # If the algorithm stops only because we have exceeded the maximum number of
        # iterations, continue improving guess of modal parameters
        total_iterations = 0
        optimization_time = 0
        max_attempts = get_setting(m, :optimization_attempts)
        attempts = 1

        while !converged
            begin_time = time_ns()
            out, H = optimize!(m, data;
                               method = get_setting(m, :optimization_method),
                               ftol = ftol, grtol = gtol, xtol = xtol,
                               iterations = n_iterations, show_trace = true, step_size = step_size,
                               verbose = verbose,
                               mle = mle)

            attempts += 1
            total_iterations += out.iterations
            converged = out.converged || attempts > max_attempts

            end_time = (time_ns() - begin_time)/1e9
            println(verbose, :low, @sprintf "Total iterations completed: %d\n" total_iterations)
            println(verbose, :low, @sprintf "Optimization time elapsed: %5.2f\n" optimization_time += end_time)

            # Write params to file after every `n_iterations` iterations
            params = map(θ->θ.value, m.parameters)
            h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
                file["params"] = params
            end
        end

        # write parameters to file one last time so we have the final mode
        h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
            file["params"] = params
        end
    end

    params = map(θ->θ.value, m.parameters)

    # Sampling does not make sense if mle=true
    if mle || !sampling
        return nothing
    end

    if get_setting(m,:sampling_method) == :MH
        ########################################################################################
        ### Step 2: Compute proposal distribution for Markov Chain Monte Carlo (MCMC)
        ###
        ### In Metropolis-Hastings, we draw sample parameter vectors from
        ### the proposal distribution, which is a degenerate multivariate
        ### normal centered at the mode. Its variance is the inverse of
        ### the hessian. We find the inverse via eigenvalue decomposition.
        ########################################################################################

        ## Calculate the Hessian at the posterior mode
        hessian = if calculate_hessian(m)
            println(verbose, :low, "Recalculating Hessian...")

            hessian, _ = hessian!(m, params, data; verbose=verbose)

            h5open(rawpath(m, "estimate","hessian.h5"),"w") do file
                file["hessian"] = hessian
            end

            hessian

    # Read in a pre-calculated Hessian
    else
        fn = hessian_path(m)
        println(verbose, :low, "Using pre-calculated Hessian from $fn")

        hessian = h5open(fn,"r") do file
            read(file, "hessian")
        end

        hessian
	end

        # Compute inverse hessian and create proposal distribution, or
        # just create it with the given cov matrix if we have it
        propdist = if isempty(proposal_covariance)
            # Make sure the mode and hessian have the same number of parameters
            n = length(params)
            @assert (n, n) == size(hessian)

            # Compute the inverse of the Hessian via eigenvalue decomposition
            S_diag, U = eigen(hessian)
            big_eig_vals = findall(x -> x > 1e-6, S_diag)
            hessian_rank = length(big_eig_vals)

            S_inv = zeros(n, n)
            for i = (n-hessian_rank+1):n
                S_inv[i, i] = 1/S_diag[i]
            end

            hessian_inv = U*sqrt.(S_inv) #this is the inverse of the hessian
            DegenerateMvNormal(params, hessian_inv)
        else
            DegenerateMvNormal(params, proposal_covariance)
        end

        if rank(propdist) != n_parameters_free(m)
            println("problem –    shutting down dimensions")
        end

        ########################################################################################
        ### Step 3: Sample from posterior using Metropolis-Hastings algorithm
        ########################################################################################

        # Set the jump size for sampling
        cc0 = get_setting(m, :mh_cc0)
        cc  = get_setting(m, :mh_cc)

        metropolis_hastings(propdist, m, data, cc0, cc; verbose = verbose,
                            filestring_addl = filestring_addl);

    elseif get_setting(m, :sampling_method) == :SMC
        ########################################################################################
        ### Step 3: Run Sequential Monte Carlo (SMC)
        ###
        ### In Sequential Monte Carlo, a large number of Markov Chains
        ### are simulated iteratively to create a particle approximation
        ### of the posterior. Portions of this method are executed in
        ### parallel.
        ########################################################################################
        smc(m, data; verbose = verbose, filestring_addl = filestring_addl,
            continue_intermediate = continue_intermediate,
            intermediate_stage_start = intermediate_stage_start,
            save_intermediate = save_intermediate,
            intermediate_stage_increment = intermediate_stage_increment)
    end

    ########################################################################################
    ### Step 4: Calculate and save parameter covariance matrix
    ########################################################################################

    compute_parameter_covariance(m, filestring_addl = filestring_addl)

    return nothing
end

"""
```
compute_parameter_covariance(m::AbstractDSGEModel)
```

Calculates the parameter covariance matrix from saved parameter draws, and writes it to the
parameter_covariance.h5 file in the `workpath(m, "estimate")` directory.

### Arguments
* `m::AbstractDSGEModel`: the model object
"""
function compute_parameter_covariance(m::AbstractDSGEModel;
                                      filestring_addl::Vector{String} = Vector{String}(undef, 0))

    sampling_method = get_setting(m, :sampling_method)
    if sampling_method ∉ [:MH, :SMC]
        throw("Invalid sampling method specified in setting :sampling_method")
    end

    prefix           = sampling_method == :MH ? "mh" : "smc"
    param_draws_path = rawpath(m, "estimate", prefix * "save.h5", filestring_addl)
    savepath         = workpath(m, "estimate", "parameter_covariance.h5", filestring_addl)

    return compute_parameter_covariance(param_draws_path, sampling_method; savepath = savepath)
end

"""
```
compute_parameter_covariance(param_draws_path::String, sampling_method::Symbol;
                             savepath = "parameter_covariance.h5")
```

Generic function calculates the parameter covariance matrix from saved parameter draws,
writes it to the specified savepath.

### Arguments
* `param_draws_path::String`: Path to file where parameter draws are stored.
* `sampling_method::Symbol`: Sampling method used for estimation.
```
   - `:MH`: Metropolis-Hastings
   - `:SMC`: Sequential Monte Carlo
```
### Optional Arguments
* `savepath::String`: Where parameter covariance matrix is to be saved. Will default
    to "parameter_covariance.h5" if unspecified.
"""
function compute_parameter_covariance(param_draws_path::String, sampling_method::Symbol;
                                      savepath::String = "parameter_covariance.h5")
    if sampling_method ∉ [:MH, :SMC]
        throw("Invalid sampling method specified in setting :sampling_method")
    end
    prefix = sampling_method == :MH ? "mh" : "smc"

    if !isfile(param_draws_path)
        @printf stderr "Saved parameter draws not found.\n"
        return
    end

    param_draws = h5open(param_draws_path, "r") do f
        read(f, prefix * "params")
    end

    # Calculate covariance matrix
    param_covariance = cov(param_draws)

    # Write to file
    h5open(savepath, "w") do f
        f[prefix * "cov"] = param_covariance
    end
end

"""
```
get_estimation_output_files(m)
```

Returns a `Dict{Symbol, String}` with all files created by `estimate(m)`.
"""
function get_estimation_output_files(m::AbstractDSGEModel)
    output_files = Dict{Symbol, String}()

    for file in [:paramsmode, :hessian, :mhsave]
        output_files[file] = rawpath(m, "estimate", "$file.h5")
    end

    for file in [:paramsmean, :parameter_covariance]
        output_files[file] = workpath(m, "estimate", "$file.h5")
    end

    for file in [:priors, :prior_posterior_means, :moments]
        output_files[file] = tablespath(m, "estimate", "$file.tex")
    end

    return output_files
end
