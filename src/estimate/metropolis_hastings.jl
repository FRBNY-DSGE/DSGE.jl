"""
```
function metropolis_hastings(propdist::Distribution,
                             loglikelihood::Function,
                             parameters::ParameterVector{S},
                             data::Matrix{T},
                             cc0::T,
                             cc::T;
                             n_blocks::Int = 1,
                             n_param_blocks::Int64 = 1,
                             n_sim::Int64          = 100,
                             n_burn::Int64         = 0,
                             mhthin::Int64         = 1,
                             adaptive_accpt::Bool  = false,
                             α::T = 1.0,      c::T = 1.0,
                             verbose::Symbol=:low,
                             savepath::String = "mhsave.h5",
                             rng::MersenneTwister = MersenneTwister(0),
                             testing::Bool = false) where {S<:Number, T<:AbstractFloat}
```

Implements the Metropolis-Hastings MCMC algorithm for sampling from the posterior
distribution of the parameters.

### Arguments

- `proposal_dist`: The proposal distribution that Metropolis-Hastings begins sampling from.
- `m`: The model object
- `data`: Data matrix for observables
- `cc0`: Jump size for initializing Metropolis-Hastings.
- `cc`: Jump size for the rest of Metropolis-Hastings.

### Optional Arguments

- `n_blocks::Int = 1`: Number of blocks (for memory-management purposes)
- `n_param_blocks::Int = 1`: Number of parameter blocks
- `n_sim::Int    = 100`: Number of simulations. Note: # saved observations will be
    n_sim * n_param_blocks * (n_blocks - b_burn).
- `n_burn::Int   = 0`: Length of burn-in period
- `mhthin::Int   = 1`: Thinning parameter (for mhthin = d, keep only every dth draw)
- `adaptive_accpt::Bool = false`: Whether or not to adaptively adjust acceptance prob.
- `α::T = 1.0`: Tuning parameter (step size) for proposal density computation in adaptive case
- `c::T = 0.5`: Tuning parameter (mixture proportion) for proposal density computation in
    adaptive case
- `verbose::Bool`: The desired frequency of function progress messages printed to
  standard out. One of:
```
   - `:none`: No status updates will be reported.
   - `:low`: Status updates provided at each block.
   - `:high`: Status updates provided at each draw.
```
- `savepath::String = "mhsave.h5"`: String specifying path to output file
- `rng::MersenneTwister = MersenneTwister(0)`: Chosen seed (overridden if testing = true)
- `testing::Bool = false`: Conditional for use when testing (determines fixed seeding)
"""
function metropolis_hastings(proposal_dist::Distribution,
                             loglikelihood::Function,
                             parameters::ParameterVector{S},
                             data::Matrix{T},
                             cc0::T,
                             cc::T;
                             n_blocks::Int64       = 1,
                             n_param_blocks::Int64 = 1,  # TODO: give these kwargs better names
                             n_sim::Int64          = 100,
                             n_burn::Int64         = 0,
                             mhthin::Int64         = 1,
                             adaptive_accpt::Bool  = false,
                             α::T = 1.0,      c::T = 0.5,
                             verbose::Symbol       = :low,
                             savepath::String      = "mhsave.h5",
                             rng::MersenneTwister  = MersenneTwister(0),
                             testing::Bool = false) where {S<:Number, T<:AbstractFloat}

    # If testing, set the random seeds at fixed numbers
    if testing
        Random.seed!(rng, 654)
    end

    propdist = init_deg_mvnormal(proposal_dist.μ, proposal_dist.σ)

    # Initialize algorithm by drawing para_old from normal distribution centered at the
    # posterior mode, until parameters within bounds (indicated by posterior value > -∞)
    para_old = rand(propdist, rng; cc = cc0)
    post_old = -Inf

    initialized = false
    while !initialized
        post_old = posterior!(loglikelihood, parameters, para_old, data; sampler = true)
        if post_old > -Inf
            propdist.μ = para_old
            initialized = true
        else
            para_old = rand(propdist, rng; cc=cc0)
        end
    end

    # Parameter Blocking
    free_para_inds = findall([!θ.fixed for θ in parameters])
    n_params       = length(parameters)
    if n_param_blocks == 1
        blocks_free = Vector{Int}[free_para_inds]
    else
        n_free_para = length(free_para_inds)
        blocks_free = SMC.generate_free_blocks(n_free_para, n_param_blocks)
        for block_f in blocks_free
            sort!(block_f)
        end
    end

    # Report number of blocks that will be used
    println(verbose, :low, "Blocks: $n_blocks")
    println(verbose, :low, "Draws per block: $n_sim")
    println(verbose, :low, "Parameter blocks: $n_param_blocks")

    # For n_sim*mhthin iterations within each block, generate a new parameter draw.
    # Decide to accept or reject, and save every (mhthin)th draw that is accepted.
    all_rejections = 0

    # Initialize matrices for parameter draws and transition matrices
    mhparams = zeros(n_sim * n_param_blocks, n_params)

    # Open HDF5 file for saving parameter draws
    simfile     = h5open(savepath, "w")
    n_saved_obs = n_sim * n_param_blocks * (n_blocks - n_burn)
    parasim     = d_create(simfile, "mhparams", datatype(Float64),
                           dataspace(n_saved_obs, n_params),
                           "chunk", (n_sim * n_param_blocks, n_params))

    # Keep track of how long metropolis_hastings has been sampling
    total_sampling_time = 0.

    for block = 1:n_blocks

        begin_time = time_ns()
        block_rejections = 0

        for j = 1:(n_sim * mhthin)

            for (k, block_a) in enumerate(blocks_free)

                # Draw para_new from the proposal distribution
                para_subset = para_old[block_a]
                d_subset    = DegenerateMvNormal(propdist.μ[block_a],
                                       (propdist.σ[block_a, block_a] +
                                       propdist.σ[block_a, block_a]') / 2.,
                                       inv((propdist.σ[block_a, block_a] +
                                       propdist.σ[block_a, block_a]') / 2.),
                                       propdist.λ_vals[block_a])

                para_draw   = rand(d_subset, rng; cc = cc)

                para_new          = deepcopy(para_old)
                para_new[block_a] = para_draw

                q0, q1 = if adaptive_accpt
                    SMC.compute_proposal_densities(para_draw, para_subset, d_subset;
                                                   α = α, c = c)
                else
                    0.0, 0.0
                end

                # Solve the model (checking that parameters are within bounds and
                # gensys returns a meaningful system) and evaluate the posterior
                post_new = posterior!(loglikelihood, parameters, para_new, data;
                                      sampler = true)

                println(verbose, :high, "Block $block, Iteration $j, Parameter Block " *
                        "$k/$(n_param_blocks): posterior = $post_new")

                # Choose to accept or reject the new parameter by calculating the
                # ratio (r) of the new posterior value relative to the old one We
                # compare min(1, r) to a number drawn randomly from a uniform (0, 1)
                # distribution. This allows us to always accept the new draw if its
                # posterior value is greater than the previous draw's, but it gives
                # some probability to accepting a draw with a smaller posterior
                # value, so that we may explore tails and other local modes.
                r = exp((post_new - post_old) + (q0 - q1))
                x = rand(rng)

                if x < min(1.0, r)
                    # Accept proposed jump
                    para_old = para_new
                    post_old = post_new
                    propdist.μ = para_new

                    println(verbose, :high, "Block $block, Iteration $j, Parameter Block " *
                        "$k/$(n_param_blocks): accept proposed jump")
                else
                    # Reject proposed jump
                    block_rejections += 1

                    println(verbose, :high, "Block $block, Iteration $j, Parameter Block " *
                        "$k/$(n_param_blocks): reject proposed jump")
                end

                # Save every (mhthin)th draw
                if j % mhthin == 0
                    draw_index = convert(Int, ((j / mhthin) - 1) * n_param_blocks + k)
                    mhparams[draw_index, :]  = para_old'
                end
            end # of loop over parameter blocks
        end # of block

        all_rejections += block_rejections
        block_rejection_rate = block_rejections / (n_sim * mhthin * n_param_blocks)

        ## Once every iblock times, write parameters to a file

        # Calculate start/end indices for this block (corresponds to new chunk in memory)
        block_start = n_sim * n_param_blocks * (block - n_burn - 1)+1
        block_end   = block_start + (n_sim * n_param_blocks) - 1

        # Write parameters to file if we're past n_burn blocks
        if block > n_burn
            parasim[block_start:block_end, :] = map(Float64, mhparams)
        end

        # Calculate time to complete this block, average block time, and
        # expected time to completion
        block_time                      = (time_ns() - begin_time) / 1e9
        total_sampling_time            += block_time
        total_sampling_time_minutes     = total_sampling_time / 60
        expected_time_remaining_sec     = (total_sampling_time / block) * (n_blocks - block)
        expected_time_remaining_minutes = expected_time_remaining_sec / 60

        println(verbose, :low, "Completed $block of $n_blocks blocks.")
        println(verbose, :low, "Total time to compute $block blocks: " *
                "$total_sampling_time_minutes minutes")
        println(verbose, :low, "Expected time remaining for Metropolis-Hastings: " *
                "$expected_time_remaining_minutes minutes")
        println(verbose, :low, "Block $block rejection rate: $block_rejection_rate \n")
    end # of loop over blocks
    close(simfile)

    rejection_rate = all_rejections / (n_blocks * n_sim * mhthin * n_param_blocks)
    println(verbose, :low, "Overall rejection rate: $rejection_rate")
end

"""
```
metropolis_hastings(propdist::Distribution, m::Union{AbstractDSGEModel,AbstractVARModel},
    data::Matrix{T}, cc0::T, cc::T; verbose::Symbol = :low) where {T<:AbstractFloat}
```

Wrapper function for DSGE models which calls Metropolis-Hastings MCMC algorithm for
sampling from the posterior distribution of the parameters.

### Arguments

- `propdist`: The proposal distribution that Metropolis-Hastings begins sampling from.
- `m`: The model object
- `data`: Data matrix for observables
- `cc0`: Jump size for initializing Metropolis-Hastings.
- `cc`: Jump size for the rest of Metropolis-Hastings.

### Optional Arguments

- `verbose`: The desired frequency of function progress messages printed to
  standard out. One of:

```
   - `:none`: No status updates will be reported.
   - `:low`: Status updates provided at each block.
   - `:high`: Status updates provided at each draw.
```
"""
function metropolis_hastings(propdist::Distribution,
                             m::Union{AbstractDSGEModel,AbstractVARModel},
                             data::Matrix{T},
                             cc0::T,
                             cc::T;
                             filestring_addl::Vector{String} = Vector{String}(undef, 0),
                             verbose::Symbol = :low) where {T<:AbstractFloat}

    n_blocks = n_mh_blocks(m)
    n_sim    = n_mh_simulations(m)
    n_burn   = n_mh_burn(m)
    mhthin   = mh_thin(m)

    n_param_blocks = n_mh_param_blocks(m)
    adaptive_accpt = get_setting(m, :mh_adaptive_accpt)
    c              = get_setting(m, :mh_c)
    α              = get_setting(m, :mh_α)

    rng      = get_rng(m)
    testing  = m.testing
    savepath = rawpath(m, "estimate", "mhsave.h5", filestring_addl)

    # To check: Defaulting to using Chandrasekhar recursions if no missing data
    use_chand_recursion = !any(isnan.(data))

    loglikelihood = if isa(m, AbstractDSGEModel)
        function _loglikelihood_dsge(p::ParameterVector, data::Matrix{Float64})::Float64
            update!(m, p)
            likelihood(m, data; sampler = true, catch_errors = false,
                       use_chand_recursion = use_chand_recursion)
        end
    elseif isa(m, AbstractVARModel)
        function _loglikelihood_var(p::ParameterVector, data::Matrix{Float64})::Float64
            update!(m, p)
            likelihood(m, data; sampler = true, catch_errors = false)
        end
    end

    return metropolis_hastings(propdist, loglikelihood, get_parameters(m), data, cc0, cc;
                               n_blocks = n_blocks, n_param_blocks = n_param_blocks,
                               adaptive_accpt = adaptive_accpt, c = c, α = α, n_sim = n_sim,
                               n_burn = n_burn, mhthin = mhthin, verbose = verbose,
                               savepath = savepath, rng = rng, testing = testing)
end
