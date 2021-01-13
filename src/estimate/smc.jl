"""
```
smc2(m::Union{AbstractDSGEModel,AbstractVARModel}, data::Matrix; verbose::Symbol, old_data::Matrix, kwargs...)
smc(m::Union{AbstractDSGEModel,AbstractVARModel},  data::DataFrame, kwargs...)
smc(m::Union{AbstractDSGEModel,AbstractVARModel},  kwargs...)
```

### Overview

These are wrapper functions to ensure simplicity of estimation of DSGE models while
navigating the DSGE package. Most keywords are passed as settings of the DSGE models
rather than keyword arguments, as in SMC.jl. For the keyword arguments that remain,
see below.

### Arguments:

- `m`: A model object, which stores parameter values, prior dists, bounds, and various
    other settings that will be referenced
- `data`: A matrix or dataframe containing the time series of the observables used in
    the calculation of the posterior/likelihood

### Keyword Arguments:
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
- `filestring_addl::Vector{String} = ""`: additional strings to be appended to the save file
    of estimation output.
- `save_intermediate::Bool = false`: Flag for whether one wants to save intermediate Cloud objects
- `intermediate_stage_increment::Int = 10`: Save Clouds at every increment
   (1 = each stage, 10 = every 10th stage, etc.). Useful if you are using a cluster with time
    limits because if you hit the time limit, then you can just
    start from an intermediate stage rather than start over.
- `continue_intermediate::Bool = false`: Flag to indicate whether one is continuing SMC from an
    intermediate stage.
- `intermediate_stage_start::Int = 0`: Intermediate stage at which one wishes to begin the estimation.
- `tempered_update_prior_weight::Float64 = 0.0,
- `run_csminwel::Bool = true`: Set to true to run the csminwel algorithm to identify the true posterior mode
    (which may not exist) after completing an estimation. The mode identified by SMC is just
    the particle with the highest posterior value, but we do not check it is actually a mode (i.e.
    the Hessian is negative definite).
- `regime_switching::Bool = false`: Set to true if there are regime-switching parameters. Otherwise, not all the values of the
    regimes will be used or saved.
- `verbose`: Desired frequency of function progress messages printed to standard out.
	- `:none`: No status updates will be reported.
	- `:low`: Status updates for SMC initialization and recursion will be included.
	- `:high`: Status updates for every iteration of SMC is output, which includes
    the mean and standard deviation of each parameter draw after each iteration,
    as well as calculated acceptance rate, ESS, and number of times resampled.

### Outputs

- `cloud`: The ParticleCloud object containing all of the information about the
    parameter values from the sample, their respective log-likelihoods, the ESS
    schedule, tempering schedule etc., which is saved in the saveroot.
"""
function smc2(m::Union{AbstractDSGEModel,AbstractVARModel}, data::Matrix{Float64};
              verbose::Symbol = :low,
              old_data::Matrix{Float64} = Matrix{Float64}(undef, size(data, 1), 0),
              old_cloud::Union{DSGE.ParticleCloud, DSGE.Cloud,
                               SMC.Cloud} = DSGE.ParticleCloud(m, 0),
              run_test::Bool = false,
              filestring_addl::Vector{String} = Vector{String}(),
              continue_intermediate::Bool = false, intermediate_stage_start::Int = 0,
              save_intermediate::Bool = false, intermediate_stage_increment::Int = 10,
              tempered_update_prior_weight::Float64 = 0.0,
              run_csminwel::Bool = true,
              regime_switching::Bool = false)

    parallel    = get_setting(m, :use_parallel_workers)
    n_parts     = get_setting(m, :n_particles)
    n_blocks    = get_setting(m, :n_smc_blocks)
    n_mh_steps  = get_setting(m, :n_mh_steps_smc)
    old_vintage = get_setting(m, :previous_data_vintage)

    smc_iteration = get_setting(m, :smc_iteration)

    λ    = get_setting(m, :λ)
    n_Φ  = get_setting(m, :n_Φ)

    # Define tempering settings
    tempered_update_prior_weight = get_setting(m, :tempered_update_prior_weight)
    tempering_target             = get_setting(m, :adaptive_tempering_target_smc)
    use_fixed_schedule           = tempering_target == 0.0

    # Step 2 (Correction) settings
    resampling_method = get_setting(m, :resampler_smc)
    threshold_ratio   = get_setting(m, :resampling_threshold)
    threshold         = threshold_ratio * n_parts

    # Step 3 (Mutation) settings
    c      = get_setting(m, :step_size_smc)
    α      = get_setting(m, :mixture_proportion)
    target = get_setting(m, :target_accept)

    use_chand_recursion = get_setting(m, :use_chand_recursion)

    my_likelihood = if isa(m, AbstractDSGEModel)
        function _my_likelihood_dsge(parameters::ParameterVector, data::Matrix{Float64})::Float64
            update!(m, parameters, regime_switching = regime_switching)
            likelihood(m, data; sampler = false, catch_errors = true,
                       use_chand_recursion = use_chand_recursion, verbose = verbose)
        end
    else isa(m, AbstractVARModel)
        function _my_likelihood_var(parameters::ParameterVector, data::Matrix{Float64})::Float64
            update!(m, parameters)
            likelihood(m, data; sampler = false, catch_errors = true, verbose = verbose)
        end
    end

    tempered_update = !isempty(old_data)

    # This step is purely for backwards compatibility purposes
    old_cloud_conv = isempty(old_cloud) ? SMC.Cloud(0,0) : SMC.Cloud(old_cloud)

    # Initialize Paths
    loadpath = ""
    if tempered_update
        if isempty(old_cloud)
            loadpath = rawpath(m, "estimate", "smc_cloud.jld2", filestring_addl)
            loadpath = replace(loadpath, r"vint=[0-9]{6}" => "vint=" * old_vintage)
        end
    elseif continue_intermediate
        loadpath = rawpath(m, "estimate", "smc_cloud", filestring_addl) *
            "_stage=$(intermediate_stage_start).jld2"
    end
    savepath = rawpath(m, "estimate", "smc_cloud.jld2", filestring_addl)
    particle_store_path = rawpath(m, "estimate", "smcsave.h5", filestring_addl)

    # Calls SMC package's generic SMC
    println("Calling SMC.jl's SMC estimation routine...")

    SMC.smc(my_likelihood, get_parameters(m), data;
            verbose = verbose,
            testing = m.testing,
            data_vintage = data_vintage(m),

            parallel = parallel,
            n_parts  = n_parts,
            n_blocks = n_blocks,
            n_mh_steps = n_mh_steps,

            λ = λ, n_Φ = n_Φ,

            resampling_method = resampling_method,
            threshold_ratio = threshold_ratio,

            c = c, α = α, target = target,

            use_fixed_schedule = use_fixed_schedule,
            tempering_target = tempering_target,

            old_data      = old_data,
            old_cloud     = old_cloud_conv,
            old_vintage   = old_vintage,
            smc_iteration = smc_iteration,

            run_test = run_test,
            filestring_addl = filestring_addl,
            loadpath = loadpath,
            savepath = savepath,
            particle_store_path = particle_store_path,

            continue_intermediate = continue_intermediate,
            intermediate_stage_start = intermediate_stage_start,
            save_intermediate = save_intermediate,
            intermediate_stage_increment = intermediate_stage_increment,
	        tempered_update_prior_weight = tempered_update_prior_weight,
            regime_switching = regime_switching)

    if run_csminwel
        m <= Setting(:sampling_method, :SMC)
        update!(m, load_draws(m, :mode, filestring_addl = filestring_addl))
        out, H = optimize!(m, data)
        println("Saving to " * replace(savepath, "smc_cloud" => "paramsmode") * "...")
        jldopen(replace(savepath, "smc_cloud" => "paramsmode"), true, true, true, IOStream) do file
            write(file, "mode", out.minimizer)
        end
    end
end

function smc(m::Union{AbstractDSGEModel,AbstractVARModel}, data::DataFrame; verbose::Symbol = :low,
             old_data::Matrix{Float64} = Matrix{Float64}(undef, size(data, 1), 0),
             old_cloud::Union{DSGE.ParticleCloud, DSGE.Cloud,
                              SMC.Cloud} = DSGE.ParticleCloud(m, 0),
             filestring_addl::Vector{String} = Vector{String}(undef, 0),
             run_test::Bool = false,
             save_intermediate::Bool = false, intermediate_stage_increment::Int = 10,
             continue_intermediate::Bool = false, intermediate_stage_start::Int = 0,
             tempered_update_prior_weight::Float64 = 0.0,
             run_csminwel::Bool = true,
             regime_switching::Bool = false)

    data_mat = df_to_matrix(m, data)
    return smc2(m, data_mat, verbose = verbose,
                old_data = old_data, old_cloud = old_cloud,
                filestring_addl = filestring_addl, run_test = run_test,
                save_intermediate = save_intermediate,
                intermediate_stage_increment = intermediate_stage_increment,
                continue_intermediate = continue_intermediate,
                intermediate_stage_start = intermediate_stage_start,
                tempered_update_prior_weight = tempered_update_prior_weight, run_csminwel = run_csminwel,
                regime_switching = regime_switching)
end

function smc(m::Union{AbstractDSGEModel,AbstractVARModel}; verbose::Symbol = :low,
             old_data::Matrix{Float64} = Matrix{Float64}(undef, size(data, 1), 0),
             old_cloud::Union{DSGE.ParticleCloud, DSGE.Cloud,
                              SMC.Cloud} = DSGE.ParticleCloud(m, 0),
             filestring_addl::Vector{String} = Vector{String}(undef, 0),
             run_test::Bool = false,
             save_intermediate::Bool = false, intermediate_stage_increment::Int = 10,
             continue_intermediate::Bool = false, intermediate_stage_start::Int = 0,
             tempered_update_prior_weight::Float64 = 0.0,
             run_csminwel::Bool = true,
             regime_switching::Bool = false)
    data = load_data(m)
    data_mat = df_to_matrix(m, data)
    return smc2(m, data_mat, verbose = verbose,
                old_data = old_data, old_cloud = old_cloud,
                filestring_addl = filestring_addl, run_test = run_test,
                save_intermediate = save_intermediate,
                intermediate_stage_increment = intermediate_stage_increment,
                continue_intermediate = continue_intermediate,
                intermediate_stage_start = intermediate_stage_start,
                tempered_update_prior_weight = tempered_update_prior_weight, run_csminwel = run_csminwel,
                regime_switching = regime_switching)
end

function isempty(c::ParticleCloud)
    length(c.particles) == 0
end

function cloud_isempty(c::ParticleCloud)
    length(c.particles) == 0
end



function get_cloud(m::AbstractDSGEModel; filepath::String = rawpath(m, "estimate", "smc_cloud.jld2"))
    return load(filepath, "cloud")
end

"""
```
function vector_particles_to_cloud(m::AbstractDSGEModel, particles::Vector{Particles})
```
Converts a ParticleCloud (old) to Cloud (new).
"""
function vector_particles_to_cloud(m::AbstractDSGEModel, particles::Vector{Particle})
    cloud = Cloud(m, length(particles))
    N = size(cloud.particles, 2)
    for i in 1:size(cloud.particles,1)
        cloud.particles[i, 1:length(m.parameters)] = particles[i].value
        cloud.particles[i, ind_loglh(N)] = particles[i].loglh
        cloud.particles[i, ind_logpost(N)] = particles[i].logpost
        update!(m, particles[i].value)
        cloud.particles[i, ind_logprior(N)] = prior(m)
        cloud.particles[i, ind_old_loglh(N)] = particles[i].old_loglh
        cloud.particles[i, ind_accept(N)] = particles[i].accept
        cloud.particles[i, ind_weight(N)] = particles[i].weight
    end
    return cloud
end
