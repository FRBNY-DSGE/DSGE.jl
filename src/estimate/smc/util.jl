"""
```
function scalar_reduce(args...)
```
Each individual iteration returns n scalars. The output is reduced to n vectors,
where the i-th vector contains all of the i-th scalars from each iteration.

The return type of reduce functions must be the same type as the tuple of
arguments passed in. If args is a tuple of Vector{Float64}, then the return
argument will be a Vector{Float64}.

e.g.
a, b = @parallel (scalar_reduce) for i in 1:10000
           [[1], [2]]
       end
a = [1, 1, 1, ...]
b = [2, 2, 2, ...]

Input/Output type: Vector{Vector{Float64}}
"""
function scalar_reduce(args...)
    return_arg = args[1]
    for (i, arg) in enumerate(args[2:end])
        for (j, el) in enumerate(arg)
            append!(return_arg[j], el)
        end
    end
    return return_arg
end

"""
```
function vector_reduce(args...)
```
Each individual iteration returns n Vector types; we vector-reduce to n matrices, where
the i-th column of that matrix corresponds to the i-th vector from an individual iteration.

Input/Output type: Vector{Matrix{Float64}}
"""
function vector_reduce(args...)
    nargs1 = length(args)    # The number of times the loop is run
    nargs2 = length(args[1]) # The number of variables output by a single run

    return_arg = args[1]
    for i in 1:nargs2
        for j in 2:nargs1
            return_arg[i] = hcat(return_arg[i], args[j][i])
        end
    end
    return return_arg
end

"""
```
function scalar_reshape(args...)
```
Function ensures type conformity of the return arguments.
"""
function scalar_reshape(args...)
    n_args = length(args)
    return_arg = Vector{Vector{Float64}}(undef, n_args)
    for i in 1:n_args
        arg = typeof(args[i]) <: Vector ? args[i] : [args[i]]
        return_arg[i] = arg
    end
    return return_arg
end

"""
```
function vector_reshape(args...)
```
Function ensures type conformity of the return arguments.
"""
function vector_reshape(args...)
    n_args = length(args)
    return_arg = Vector{Matrix{Float64}}(undef, n_args)
    for i in 1:n_args
        arg = typeof(args[i]) <: Vector ? args[i] : [args[i]]
        return_arg[i] = reshape(arg, length(arg), 1)
    end
    return return_arg
end

"""
```
function generate_free_blocks(n_free_para, n_blocks)
```

Return a Vector{Vector{Int64}} where each internal Vector{Int64} contains a subset of the range 1:n_free_para of randomly permuted indices. This is used to index out random blocks of free parameters from the covariance matrix for the mutation step.
"""
function generate_free_blocks(n_free_para::Int64, n_blocks::Int64)
    rand_inds = shuffle(1:n_free_para)

    subset_length     = cld(n_free_para, n_blocks) # ceiling division
    last_block_length = n_free_para - subset_length*(n_blocks - 1)

    blocks_free = Vector{Vector{Int64}}(undef, n_blocks)
    for i in 1:n_blocks
        if i < n_blocks
            blocks_free[i] = rand_inds[((i-1)*subset_length + 1):(i*subset_length)]
        else
            # To account for the fact that the last block may be smaller than the others
            blocks_free[i] = rand_inds[end-last_block_length+1:end]
        end
    end
    return blocks_free
end

function isempty(c::ParticleCloud)
    length(c.particles) == 0
end

"""
```
function generate_all_blocks(blocks_free, free_para_inds)
```

Return a Vector{Vector{Int64}} where each internal Vector{Int64} contains indices corresponding to those in `blocks_free` but mapping to `1:n_para` (as opposed to `1:n_free_para`). These blocks are used to reconstruct the particle vector by inserting the mutated free parameters into the size `n_para,` particle vector, which also contains fixed parameters.
"""
function generate_all_blocks(blocks_free::Vector{Vector{Int64}}, free_para_inds::Vector{Int64})
    n_free_para = length(free_para_inds)
    ind_mappings = Dict{Int64, Int64}()

    for (k, v) in zip(1:n_free_para, free_para_inds)
        ind_mappings[k] = v
    end

    blocks_all = similar(blocks_free)
    for (i, block) in enumerate(blocks_free)
        blocks_all[i] = similar(block)
        for (j, b) in enumerate(block)
            blocks_all[i][j] = ind_mappings[b]
        end
    end
    return blocks_all
end

function get_cloud(m::AbstractModel; filepath::String = rawpath(m, "estimate", "smc_cloud.jld2"))
    return load(filepath, "cloud")
end

function init_stage_print(cloud::ParticleCloud;
                          verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    if use_fixed_schedule
        println("--------------------------")
            println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
    else
        println("--------------------------")
            println("Iteration = $(cloud.stage_index)")
    end
	println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
	println("--------------------------")
        println("c = $(cloud.c)")
        println("ESS = $(cloud.ESS[cloud.stage_index])   ($(cloud.resamples) total resamples.)")
	println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], digits = 5)), $(round(σ[n], digits = 5))")
	    end
    end
end

function init_stage_print(cloud::Cloud, para_symbols::Vector{Symbol};
                          verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    if use_fixed_schedule
        println("--------------------------")
            println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
    else
        println("--------------------------")
            println("Iteration = $(cloud.stage_index)")
    end
	println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
	println("--------------------------")
        println("c = $(cloud.c)")
        println("ESS = $(cloud.ESS[cloud.stage_index])   ($(cloud.resamples) total resamples.)")
	println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(para_symbols)
            println("$(para_symbols[n]) = $(round(μ[n], digits = 5)), $(round(σ[n], digits = 5))")
	    end
    end
end

function end_stage_print(cloud::ParticleCloud;
                         verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    total_sampling_time_minutes = cloud.total_sampling_time/60
    if use_fixed_schedule
        expected_time_remaining_sec = (cloud.total_sampling_time/cloud.stage_index)*(cloud.n_Φ - cloud.stage_index)
        expected_time_remaining_minutes = expected_time_remaining_sec/60
    end

    println("--------------------------")
    if use_fixed_schedule
        println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
        println("time elapsed: $(round(total_sampling_time_minutes, digits = 4)) minutes")
        println("estimated time remaining: $(round(expected_time_remaining_minutes, digits = 4)) minutes")
    else
        println("Iteration = $(cloud.stage_index)")
        println("time elapsed: $(round(total_sampling_time_minutes, digits = 4)) minutes")
    end
    println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
    println("--------------------------")
        println("c = $(cloud.c)")
        println("accept = $(cloud.accept)")
        println("ESS = $(cloud.ESS[cloud.stage_index])   ($(cloud.resamples) total resamples.)")
    println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(cloud.particles[1])
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], digits = 5)), $(round(σ[n], digits = 5))")
        end
    end
end
function end_stage_print(cloud::Cloud, para_symbols::Vector{Symbol};
                         verbose::Symbol=:low, use_fixed_schedule::Bool = true)
    total_sampling_time_minutes = cloud.total_sampling_time/60
    if use_fixed_schedule
        expected_time_remaining_sec = (cloud.total_sampling_time/cloud.stage_index)*(cloud.n_Φ - cloud.stage_index)
        expected_time_remaining_minutes = expected_time_remaining_sec/60
    end

    println("--------------------------")
    if use_fixed_schedule
        println("Iteration = $(cloud.stage_index) / $(cloud.n_Φ)")
        println("time elapsed: $(round(total_sampling_time_minutes, digits = 4)) minutes")
        println("estimated time remaining: $(round(expected_time_remaining_minutes, digits = 4)) minutes")
    else
        println("Iteration = $(cloud.stage_index)")
        println("time elapsed: $(round(total_sampling_time_minutes, digits = 4)) minutes")
    end
    println("--------------------------")
        println("phi = $(cloud.tempering_schedule[cloud.stage_index])")
    println("--------------------------")
        println("c = $(cloud.c)")
        println("accept = $(cloud.accept)")
        println("ESS = $(cloud.ESS[cloud.stage_index])   ($(cloud.resamples) total resamples.)")
    println("--------------------------")
    if VERBOSITY[verbose] >= VERBOSITY[:high]
        μ = weighted_mean(cloud)
        σ = weighted_std(cloud)
        for n=1:length(para_symbols)
            println("$(para_symbols[n]) = $(round(μ[n], digits = 5)), $(round(σ[n], digits = 5))")
        end
    end
end
