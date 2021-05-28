"""
```
estimate(m, data; verbose=:low, proposal_covariance=Matrix())
```

Estimate the DSGE parameter posterior distribution.

### Arguments:
- `m::AbstractDSGEModel` or `m::AbstractVARModel`: model object

### Estimation Settings
Please see the section on 'Estimation Settings' on the 'Advanced Usage'
page of the online documentation
or src/defaults.jl for a full description of
all the estimation settings for both
the Metropolis-Hastings (MH) and Sequential Monte Carlo (SMC) algorithms.
For the latter, also see `?DSGE.smc2`. Most of the
optional and keyword arguments described below are
not directly related to the behavior of the sampling algorithms
(e.g. tuning, number of samples).

### Optional Arguments:
- `data`: well-formed data as `Matrix` or `DataFrame`. If this is not provided, the `load_data` routine will be executed.

### Keyword Arguments:
- `verbose::Symbol`: The desired frequency of function progress messages printed to standard out.
   - `:none`: No status updates will be reported.
   - `:low` (default): Status updates will be provided in csminwel and at each block in
     Metropolis-Hastings. For SMC's verbosity settings, see `?smc2`.
   - `:high`: Status updates provided at each iteration in Metropolis-Hastings.
- `proposal_covariance::Matrix = []`: Used to test the metropolis_hastings algorithm with a precomputed
  covariance matrix for the proposal distribution. When the Hessian is singular,
  eigenvectors corresponding to zero eigenvectors are not well defined, so eigenvalue
  decomposition can cause problems. Passing a precomputed matrix allows us to ensure that
  the rest of the routine has not broken.
- `method::Symbol`: The method to use when sampling from the posterior distribution. Can
    be either `:MH` for standard Metropolis Hastings Markov Chain Monte Carlo, or `:SMC`
    for Sequential Monte Carlo. This should be specified by the setting `sampling_method` in `m`.
- `mle = false`: Set to true if parameters should be estimated by maximum likelihood directly.
    If this is set to true, this function will return after estimating parameters.
- `sampling = true`: Set to false to disable sampling from the posterior.
- `old_data::Matrix{Float64}` = []: A matrix containing the time series of observables of previous data
    (with `data` being the new data) for the purposes of a time tempered estimation
    (that is, using the posterior draws from a previous estimation as the initial set
    of draws for an estimation with new data). Running a bridge estimation
    requires both `old_data` and `old_cloud`.
- `old_cloud::Union{DSGE.ParticleCloud, DSGE.Cloud, SMC.Cloud} = DSGE.Cloud(m, 0)`: old Cloud object
    used to describe particles from a previous estimation with old data. Running a bridge estimation
    requires both `old_data` and `old_cloud`. If running a bridge estimation and
    no `old_cloud` is provided, then it will be loaded using
    the filepaths in `m`. If no `old_cloud` exists, then the bridge estimation will not run.
- `old_model::Union{AbstractDSGEModel, AbstractVARModel} = m`: model object from which we can build
    the old log-likelihood function for a time tempered SMC estimation. It should be possible
    to evaluate the old log-likelihood given `old_data` and the current draw of parameters.
    This may be nontrivial if, for example, *new* parameters have been added to `m` since the old
    estimation. In this case, `old_model` should include the *new* parameters but still return
    the old log-likelihood as the original estimation if given the same *old* parameters.
    By default, we assume the log-likelihood function has not changed
    and therefore coincides with the current one.
- `filestring_addl::Vector{String} = []`: Additional strings to add to the file name
    of estimation output as a way to distinguish output from each other.
- `continue_intermediate::Bool = false`: set to true if the estimation is starting
    from an intermediate stage that has been previously saved.
- `intermediate_stage_start::Int = 0`: number of the stage from which the user wants
    to continue the estimation (see `continue_intermediate`)
- `save_intermediate::Bool = true`: set to true to save intermediate stages when using SMC
- `intermediate_stage_increment::Int = 10`: number of stages that must pass before saving
    another intermediate stage.
-  `run_csminwel::Bool = true`: by default, csminwel is run after a SMC estimation finishes
    to recover the true mode of the posterior. Set to false to avoid this step
    (csminwel can take hours for medium-scale DSGE models).
-  `toggle::Bool = true`: when regime-switching, several functions assume regimes
    are toggled to regime 1. If the likelihood function is not
    written to toggle to regime 1 when done, then regime-switching estimation
    will not work properly. Set to `false` to reduce computation time if the
    user is certain that the likelihood is written properly.
"""
function estimate(m::Union{AbstractDSGEModel,AbstractVARModel}, df::DataFrame;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0, 0),
                  mle::Bool = false,
                  sampling::Bool = true,
                  filestring_addl::Vector{String} = Vector{String}(),
                  old_data::Matrix{Float64} = Matrix{Float64}(undef, size(data, 1), 0),
                  old_cloud::Union{DSGE.ParticleCloud, DSGE.Cloud,
                                   SMC.Cloud} = DSGE.ParticleCloud(m, 0),
                  old_model::Union{AbstractDSGEModel, AbstractVARModel} = m,
                  continue_intermediate::Bool = false,
                  intermediate_stage_start::Int = 0,
                  intermediate_stage_increment::Int = 10,
                  save_intermediate::Bool = false,
                  run_csminwel::Bool = true,
                  toggle::Bool = true,
                  log_prob_oldy::Float64 = 0.0)
    data = df_to_matrix(m, df)
    estimate(m, data; verbose = verbose, proposal_covariance = proposal_covariance,
             mle = mle, sampling = sampling,
             old_data = old_data, old_cloud = old_cloud,
             old_model = old_model,
             continue_intermediate = continue_intermediate,
             intermediate_stage_increment = intermediate_stage_increment,
             save_intermediate = save_intermediate,
             run_csminwel = run_csminwel, toggle = toggle, log_prob_oldy = log_prob_oldy)
end

function estimate(m::Union{AbstractDSGEModel,AbstractVARModel};
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0, 0),
                  mle::Bool = false,
                  sampling::Bool = true,
                  filestring_addl::Vector{String} = Vector{String}(),
                  old_data::Matrix{Float64} = Matrix{Float64}(undef, size(data, 1), 0),
                  old_cloud::Union{DSGE.ParticleCloud, DSGE.Cloud,
                                   SMC.Cloud} = DSGE.ParticleCloud(m, 0),
                  old_model::Union{AbstractDSGEModel, AbstractVARModel} = m,
                  continue_intermediate::Bool = false,
                  intermediate_stage_start::Int = 0,
                  intermediate_stage_increment::Int = 10,
		          save_intermediate::Bool = false,
                  run_csminwel::Bool = true,
                  toggle::Bool = true, log_prob_oldy::Float64 = 0.0)
    # Load data
    df = load_data(m; verbose = verbose)
    estimate(m, df; verbose = verbose, proposal_covariance = proposal_covariance,
             mle = mle, sampling = sampling,
             old_data = old_data, old_cloud = old_cloud,
             old_model = old_model,
             continue_intermediate = continue_intermediate,
             intermediate_stage_increment = intermediate_stage_increment,
	         save_intermediate = save_intermediate,
             run_csminwel = run_csminwel, toggle = toggle, log_prob_oldy = log_prob_oldy)
end

function estimate(m::Union{AbstractDSGEModel,AbstractVARModel}, data::AbstractArray;
                  verbose::Symbol = :low,
                  proposal_covariance::Matrix = Matrix(undef, 0,0),
                  mle::Bool = false,
                  sampling::Bool = true,
                  filestring_addl::Vector{String} = Vector{String}(),
                  old_data::Matrix{Float64} = Matrix{Float64}(undef, size(data, 1), 0),
                  old_cloud::Union{DSGE.ParticleCloud, DSGE.Cloud,
                                   SMC.Cloud} = DSGE.ParticleCloud(m, 0),
                  old_model::Union{AbstractDSGEModel, AbstractVARModel} = m,
                  continue_intermediate::Bool = false,
                  intermediate_stage_start::Int = 0,
                  intermediate_stage_increment::Int = 10,
		          save_intermediate::Bool = false,
                  run_csminwel::Bool = true,
                  toggle::Bool = true, log_prob_oldy::Float64 = 0.0)

    if !(get_setting(m, :sampling_method) in [:SMC, :MH])
        error("method must be :SMC or :MH")
    else
        method = get_setting(m, :sampling_method)
    end

    regime_switching = haskey(get_settings(m), :regime_switching) &&
        get_setting(m, :regime_switching)

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
                               mle = mle, toggle = toggle, verbose = verbose)


            attempts += 1
            total_iterations += out.iterations
            converged = out.converged || attempts > max_attempts

            end_time = (time_ns() - begin_time)/1e9
            println(verbose, :low, @sprintf "Total iterations completed: %d\n" total_iterations)
            println(verbose, :low, @sprintf "Optimization time elapsed: %5.2f\n" optimization_time += end_time)

            # Write params to file after every `n_iterations` iterations
            params = ModelConstructors.get_values(get_parameters(m); regime_switching = regime_switching)
            h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
                file["params"] = params
            end
        end

        # write parameters to file one last time so we have the final mode
        h5open(rawpath(m, "estimate", "paramsmode.h5"),"w") do file
            file["params"] = params
        end
    end

    params = ModelConstructors.get_values(get_parameters(m); regime_switching = regime_switching)

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

            hessian, _ = hessian!(m, params, data; toggle = toggle, verbose = verbose)

            h5open(rawpath(m, "estimate","hessian.h5"),"w") do file
                file["hessian"] = hessian
            end

            hessian

        ## Read in a pre-calculated Hessian
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
            #S_diag, U = eigen(hessian)
            F = svd(hessian)

            big_eig_vals = findall(x -> x > 1e-6, F.S)
            hessian_rank = length(big_eig_vals)

            S_inv = zeros(n, n)
            #for i = (n-hessian_rank+1):n
            for i = 1:hessian_rank
                S_inv[i, i] = 1/F.S[i]
            end

            #hessian_inv = U*sqrt.(S_inv) # this is the inverse of the hessian
            hessian_inv = F.V * S_inv * F.U'#sqrt.(S_inv) * F.U'
            DegenerateMvNormal(params, hessian_inv, hessian, diag(S_inv))
        else
            DegenerateMvNormal(params, proposal_covariance, pinv(proposal_covariance),
                               eigen(proposal_covariance).values)
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

        metropolis_hastings(propdist, m, data, cc0, cc; regime_switching = regime_switching,
                            toggle = toggle, verbose = verbose, filestring_addl = filestring_addl);

    elseif get_setting(m, :sampling_method) == :SMC
        ########################################################################################
        ### Step 3: Run Sequential Monte Carlo (SMC)
        ###
        ### In Sequential Monte Carlo, a large number of Markov Chains
        ### are simulated iteratively to create a particle approximation
        ### of the posterior. Portions of this method are executed in
        ### parallel.
        ########################################################################################
        smc2(m, data; verbose = verbose, filestring_addl = filestring_addl,
             old_data = old_data, old_cloud = old_cloud,
             old_model = old_model,
             continue_intermediate = continue_intermediate,
             intermediate_stage_start = intermediate_stage_start,
             save_intermediate = save_intermediate,
             intermediate_stage_increment = intermediate_stage_increment,
             run_csminwel = run_csminwel,
             regime_switching = regime_switching, log_prob_oldy = log_prob_oldy)
    end

    ########################################################################################
    ### Step 4: Calculate and save parameter covariance matrix
    ########################################################################################

    compute_parameter_covariance(m, filestring_addl = filestring_addl)

    return nothing
end

"""
```
compute_parameter_covariance(m::Union{AbstractDSGEModel,AbstractVARModel})
```

Calculates the parameter covariance matrix from saved parameter draws, and writes it to the
parameter_covariance.h5 file in the `workpath(m, "estimate")` directory.

### Arguments
* `m::Union{AbstractDSGEModel,AbstractVARModel}`: the model object
"""
function compute_parameter_covariance(m::Union{AbstractDSGEModel,AbstractVARModel};
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
function get_estimation_output_files(m::Union{AbstractDSGEModel,AbstractVARModel})
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
