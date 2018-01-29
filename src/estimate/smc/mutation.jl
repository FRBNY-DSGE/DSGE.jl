"""
```
mutation(m, data, p, std_dev_mat, ϕ_n, ϕ_n1; c = 1., p_prop)
```

Execute one proposed move of the Metropolis-Hastings algorithm for a given parameter

### Arguments:
- `m::AbstractModel`: Model of type AbstractModel being estimated.
- `data::Matrix{Float64}`: Matrix of data
- `p::Particle`: Initial particle value
- `std_dev_mat::Matrix{Float64}`: The std_dev matrix of all initial particle values (taken w/ SVD)
- `ϕ_n::Float64`: The current tempering factor
- `ϕ_n1::Float64`: The previous tempering factor

### Keyword Arguments:
- `c::Float64`: The scaling parameter for the proposed covariance matrix
- `p_prop::Vector{Float64}`: The proposed particle value to be mixed in the mvnormal_mixture_draw
- `α::Float64`: The mixing proportion
- `old_data::Matrix{Float64}`: The matrix of old data to be used in calculating the old_loglh, old_logpost in time tempering

### Outputs:

- `p::Particle`: An updated particle containing updated parameter values, log-likelihood, posterior, and acceptance indicator.

"""
function mutation(m::AbstractModel, data::Matrix{Float64}, p::Particle, std_dev_mat::Matrix{Float64},
                  ϕ_n::Float64, ϕ_n1::Float64; c::Float64 = 1., p_prop::Vector{Float64} = p.value,
                  α::Float64 = 1., old_data::Matrix{Float64} = Matrix{Float64}(size(data, 1), 0))

    n_para = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_steps = get_setting(m, :n_MH_steps_smc)

    fixed_para_inds = find([θ.fixed for θ in m.parameters])
    nonfixed_para_inds = find([!θ.fixed for θ in m.parameters])

    # draw initial step probability
    # conditions for testing purposes
    step_prob = rand()

    # Ensuring numerical error does not propagate and move fixed parameters
    std_dev_mat[:, fixed_para_inds] = std_dev_mat[fixed_para_inds, :] = 0.

    para = p.value
    like = p.loglh
    post_init = p.logpost
    like_prev = p.old_loglh # The likelihood evaluated at the old data (for time tempering)

    # Previous posterior needs to be updated (due to tempering)
    post = post_init + (ϕ_n - ϕ_n1) * like
    accept = false

    ### BLOCKING ###

    # Generate random ordering
    # Returns indices sorted by their random value
    para_inds = [rand() for j in 1:n_para]
    para_inds = sortperm(para_inds)

    # Find the indices for the lower and upper bounds on each block
    b = 1
    blo = Array{Int}(n_blocks)
    bhi = Array{Int}(n_blocks)
    for b in 1:n_blocks
        blo[b] = convert(Int, max(0, ceil((b-1)*n_para/n_blocks)))
        bhi[b] = convert(Int, min(n_para+1, ceil(b*n_para/n_blocks)+1))
    end

    l = 0
    while l < n_steps
        for (j, k) in zip(blo, bhi)
            #Zeroing out relevant rows in cov matrix
            std_dev_mat_block = copy(std_dev_mat)
            blo_inds = para_inds[1:j]
            bhi_inds = para_inds[k:end]

            std_dev_mat_block[:,blo_inds] = 0.
            std_dev_mat_block[blo_inds,:] = 0.
            std_dev_mat_block[:,bhi_inds] = 0.
            std_dev_mat_block[bhi_inds,:] = 0.

            para_new, para_new_density, para_old_density = mvnormal_mixture_draw(para, std_dev_mat_block;
                                                                                 cc = c, α = α,
                                                                                 θ_prop = p_prop)
            like_new = -Inf
            post_new = -Inf
            like_old_data = -Inf
            try
                update!(m, para_new)
                like_new = likelihood(m, data; sampler = true)
                post_new = ϕ_n*like_new + prior(m) - para_new_density
                like_old_data = likelihood(m, old_data; sampler = true)
            catch
                post_new = like_new = like_old_data = -Inf
            end

            # Accept/Reject
            post_old = post - para_old_density
            η = exp(post_new - post_old)

            if step_prob < η # accept
                para   = para_new
                like   = like_new
                post   = post_new + para_new_density # Have to add it back so as to not accumulate
                like_prev = like_old_data            # para_new_density throughout the iterations
                accept = true
            end

            # draw again for the next step
            step_prob = rand()
        end
        l += 1
    end
    update_mutation!(p, para, like, post, like_prev, accept)
    return p
end
