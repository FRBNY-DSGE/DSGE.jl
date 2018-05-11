# For mutation
# Generate a Vector of Vector{Int64} of length n_blocks, where each
# element contains a subset of the randomly permuted set of indices 1:n_para
# For the purpose of indexing into R_fr
function generate_free_blocks(n_para::Int64, n_blocks::Int64)
    rand_inds = shuffle(1:n_para)

    subset_length = cld(n_para, n_blocks) # ceiling division
    last_block_length = n_para - subset_length*(n_blocks - 1)

    blocks_free = Vector{Vector{Int64}}(n_blocks)
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

# For mutation
# Generate a Vector of Vector{Int64} of length n_blocks, where each
# element contains a subset corresponding to the subset of blocks_free of the same
# index but with indices that map to free_para_inds as opposed to 1:n_para
# For the purpose of "re-creating" the proposed parameter vector that contains both free
# and fixed parameters from the mh step generated from only the free parameters
function generate_all_blocks(blocks_free::Vector{Vector{Int64}}, free_para_inds::Vector{Int64})
    n_free_para = length(free_para_inds)
    # Need to know the mapping from an ordered list of 1:n_free_para
    # to the index in the actual parameter vector
    ind_mappings = Dict{Int64, Int64}()
    for (k, v) in zip(1:n_free_para, free_para_inds)
        ind_mappings[k] = v
    end

    # Want: Input: blocks, a vector of vectors of indices of randomized blocks of an ordered list of 1:n_free_para
    # Output: rev_blocks, a vector of vector of indices of the rand blocks indices that correspond to the actual
    # indices in a parameter vector
    function block_map(blocks::Vector{Vector{Int64}}, ind_mappings::Dict{Int64, Int64})
        blocks_all = similar(blocks)
        for (i, block) in enumerate(blocks)
            blocks_all[i] = similar(block)
            for (j, b) in enumerate(block)
                blocks_all[i][j] = ind_mappings[b]
            end
        end
        return blocks_all
    end

    blocks_all = block_map(blocks_free, ind_mappings)

    return blocks_all
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
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
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
        println("time elapsed: $(round(total_sampling_time_minutes, 4)) minutes")
        println("estimated time remaining: $(round(expected_time_remaining_minutes, 4)) minutes")
    else
        println("Iteration = $(cloud.stage_index)")
        println("time elapsed: $(round(total_sampling_time_minutes, 4)) minutes")
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
            println("$(cloud.particles[1].keys[n]) = $(round(μ[n], 5)), $(round(σ[n], 5))")
        end
    end
end
