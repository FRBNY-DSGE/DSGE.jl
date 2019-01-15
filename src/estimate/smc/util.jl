# Utility function for checking all of the files in a directory. Want to grab when resamples occur.
function get_resample_periods(path::String)
    resample_periods = Base.filter(x -> (endswith(x,"resamp.txt") &
                                         (x != "resamp.txt")), readdir(path))
    return [parse(Int, file[1:3]) for file in resample_periods]
end

function resample(period::Int, path::String)
    return Int.(vec(readdlm(path * convert_string(period) * "resamp.txt")))
end

function convert_string(i::Int)
    if i < 10
        return "00" * string(i)
    elseif i < 100
        return "0" * string(i)
    else
        return string(i)
    end
end

function parameter_map()
    return Dict{Int64, Int64}(
        1 => 7,   # varphi => :S''
        2 => 19,  # sigmac => :σ_c
        3 => 8,   # h      => :h
        4 => 13,  # lamw   => :λ_w    || FIXED
        5 => 11,  # xiw    => :ζ_w
        6 => 22,  # epsw   => :ϵ_w    || FIXED
        7 => 10,  # sigmal => :ν_l
        8 => 4,   # del    => :δ      || FIXED
        9 => 2,   # xip    => :ζ_p
        10 => 21, # epsp   => :ϵ_p    || FIXED
        11 => 12, # iotaw  => :ι_w
        12 => 3,  # iotap  => :ι_p
        13 => 9,  # ppsi   => :ppsi
        14 => 6,  # capphi => :Φ
        15 => 15, # rpi    => :ψ1
        16 => 20, # rho    => :ρ
        17 => 16, # ry     => :ψ2
        18 => 17, # rdely  => :ψ3
        19 => 18, # pibar  => :π_star
        20 => 14, # betin  => :β
        21 => 24, # lbar   => :Lmean   // given different prior
        22 => 23, # gambar => :γ
        23 => 1,  # alp    => :α
        24 => 25, # g_y    => :g_star  || FIXED

        25 => 29, # rhoa   => :ρ_z
        26 => 27, # rhob   => :ρ_b
        27 => 26, # rhog   => :ρ_g
        28 => 28, # rhoi   => :ρ_μ
        29 => 32, # rhor   => :ρ_rm
        30 => 30, # rhop   => :ρ_λ_f
        31 => 31, # rhow   => :ρ_λ_w

        32 => 41, # mup    => :η_λ_f
        33 => 42, # muw    => :η_λ_w
        34 => 40, # rhoga  => :η_gz

        35 => 36, # sigmaa => :σ_z
        36 => 34, # sigmab => :σ_b
        37 => 33, # sigmag => :σ_g
        38 => 35, # sigmai => :σ_μ
        39 => 39, # sigmar => :σ_rm
        40 => 37, # sigmap => :σ_λ_f
        41 => 38) # sigmaw => :σ_λ_w
end

# The return type of reduce functions must be the same type as the tuple of arguments being input
# E.g. If args is a tuple of Vector{Float64}, then the return argument must also be a Vector{Float64}
# Thus, to implement a scalar reduce function, where each individual iteration returns
# n scalars, and we want the output to be reduced to n vectors, where the i-th vector
# contains all of the i-th scalars from each individual iteration, then we must modify the
# individual iterations to return n singleton vectors (one element vectors) of those n
# scalars so as to preserve the homogeneity of the input/output type coming into/out of
# the scalar reduce function.
# This would not work if the input argument types were just Ints or Float64s
# since a return type of Int/Float64 for scalar reduce function does not permit that
# function to return collections of items (because an Int/Float64 can only contain a
# single value).
# E.g.
# a, b = @distributed (scalar_reduce) for i in 1:10000
    # [[1], [2]]
# end
# a = [1, 1, 1, ...]
# b = [2, 2, 2, ...]

# Input/Output type: Vector{Vector{Float64}}
function scalar_reduce(args...)
    return_arg = args[1]
    for (i, arg) in enumerate(args[2:end])
        for (j, el) in enumerate(arg)
            append!(return_arg[j], el)
        end
    end
    return return_arg
end

# Same logic applies to the vector reduce, where each individual iteration returns n
# Vector types, and we want to vector reduce to n matrices, where the i-th column of that
# matrix corresponds to the i-th vector from an individual iteration.
# Input/Output type: Vector{Matrix{Float64}}
function vector_reduce(args...)
    nargs1 = length(args) # The number of times the loop is run
    nargs2 = length(args[1]) # The number of variables output by a single run

    return_arg = args[1]
    for i in 1:nargs2
        for j in 2:nargs1
            return_arg[i] = hcat(return_arg[i], args[j][i])
        end
    end
    return return_arg
end

# The following two functions are to ensure type conformity of the return arguments
function scalar_reshape(args...)
    n_args = length(args)
    return_arg = Vector{Vector{Float64}}(undef, n_args)
    for i in 1:n_args
        arg = typeof(args[i]) <: Vector ? args[i] : [args[i]]
        return_arg[i] = arg
    end
    return return_arg
end

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

    subset_length = cld(n_free_para, n_blocks) # ceiling division
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

"""
```
function generate_all_blocks(n_free_para, n_blocks, path, i)
```
USED FOR FORTRAN TESTING.
Return a Vector{Vector{Int64}} where each internal Vector{Int64} contains a subset of the range 1:n_free_para of randomly permuted indices. This is used to index out random blocks of free parameters from the covariance matrix for the mutation step.
"""
function generate_all_blocks(n_free_para::Int64, n_blocks::Int64, path::String, i::Int64)
    if n_blocks == 3
        break_points = [0, 13, 25, 36] # this is from the FORTRAN code
    else #n_blocks == 1
        break_points = [0, 36]
    end
        # fort_inds    = readdlm(path * convert_string(i) * "randomblocks.txt")[1:n_free_para]
    # my_map       = parameter_map()
    # inds         = [my_map[x] for x in fort_inds]
    inds = readdlm(path * convert_string(i) * "randomblocks.txt")[1:n_free_para]
    inds = [Int(x) for x in inds]
    fortran_list = [inds[break_points[j]+1:break_points[j+1]] for j=1:n_blocks]
    return fortran_list
end

"""
```
function generate_free_blocks(blocks_free, free_para_inds)
```
FOR FOTRAN TESTING.
Return a Vector{Vector{Int64}} where each internal Vector{Int64} contains indices corresponding to those in `blocks_free` but mapping to `1:n_para` (as opposed to `1:n_free_para`). These blocks are used to reconstruct the particle vector by inserting the mutated free parameters into the size `n_para,` particle vector, which also contains fixed parameters.
"""
function generate_free_blocks(blocks_all::Vector{Vector{Int64}}, free_para_inds::Vector{Int64})
    blocks_free = similar(blocks_all)
    my_map = Dict{Int64, Int64}()
    for (free_ind, ind) in zip([i for i=1:length(free_para_inds)], free_para_inds)
        my_map[ind] = free_ind
    end
    for (i, block) in enumerate(blocks_all)
        blocks_free[i] = [my_map[x] for x in block]
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

function get_cloud(m::AbstractModel; filepath::String = rawpath(m, "estimate", "smc_cloud.jld"))
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
