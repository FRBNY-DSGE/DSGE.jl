"""
```
mutation(m, data, p, R, ϕ_n, ϕ_n1; c = 1., p_prop)
```

Execute one proposed move of the Metropolis-Hastings algorithm for a given parameter

### Arguments:
- `m::AbstractModel`: Model of type AbstractModel being estimated.
- `data::Matrix{Float64}`: Matrix of data
- `p::Particle`: Initial particle value
- `R::Matrix{Float64}`: The covariance matrix of all initial particle values
- `ϕ_n::Float64`: The current tempering factor
- `ϕ_n1::Float64`: The previous tempering factor

### Keyword Arguments:
- `c::Float64`: The scaling parameter for the proposed covariance matrix
- `p_prop::Vector{Float64}`: The proposed particle value to be mixed in the mvnormal_mixture_draw

### Outputs:

- `p::Particle`: An updated particle containing updated parameter values, log-likelihood, posterior, and acceptance indicator.

"""
function mutation(m::AbstractModel, data::Matrix{Float64}, p::Particle, R::Matrix{Float64},
                  ϕ_n::Float64, ϕ_n1::Float64; c::Float64 = 1., p_prop::Vector{Float64} = p.value)

    n_para = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_steps = get_setting(m, :n_MH_steps_smc)

    fixed_para_inds = find([ θ.fixed for θ in m.parameters])

    # draw initial step probability
    # conditions for testing purposes
    step_prob = rand()

    # This calculates the matrix 'square root' of R for the purposes
    # of drawing by μ + std_dev_mat*randn(length(μ))
    U, E, V = svd(R)
    std_dev_mat = U * diagm(sqrt.(E))

    # Ensuring numerical error does not propagate and move fixed parameters
    std_dev_mat[:, fixed_para_inds] = std_dev_mat[fixed_para_inds, :] = 0.

    para = p.value
    like = p.loglh
    post_init = p.logpost

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

            std_dev_mat_block[:,blo_inds] = 0
            std_dev_mat_block[blo_inds,:] = 0
            std_dev_mat_block[:,bhi_inds] = 0
            std_dev_mat_block[bhi_inds,:] = 0

            para_new, para_new_density, para_old_density = mvnormal_mixture_draw(para, std_dev_mat_block;
                                                                                 cc = c, α = .9,
                                                                                 θ_prop = p_prop)
            like_new = -Inf
            post_new = -Inf
            try
                update!(m, para_new)
                like_new = likelihood(m, data; sampler = true)
                post_new = ϕ_n*like_new + prior(m) - para_new_density
            catch err
                if isa(err, GensysError) || isa(err, ParamBoundsError)
                    post_new = like_new = -Inf
                # Ignoring cases where G0 in gensys has a unit root
                elseif isa(err, Base.LinAlg.LAPACKException)
                    post_new = like_new = -Inf
                else
                    throw(err)
                end
            end

            # Accept/Reject
            post_old = post - para_old_density
            α = exp(post_new - post_old)

            if step_prob < α # accept
                para   = para_new
                like   = like_new
                post   = post_new + para_new_density # Have to add it back so as to not accumulate 
                accept = true                        # para_new_density throughout the iterations
            end

            # draw again for the next step
            step_prob = rand()
        end
        l += 1
    end
    update_mutation!(p, para, like, post, accept)
    return p
end
