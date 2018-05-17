"""
```
mutation(m, data, p, d, blocks_free, blocks_all, ϕ_n, ϕ_n1; c = 1., α = 1., old_data)
```

Execute one proposed move of the Metropolis-Hastings algorithm for a given parameter

### Arguments:
- `m::AbstractModel`: Model of type AbstractModel being estimated.
- `data::Matrix{Float64}`: Matrix of data
- `p::Particle`: Initial particle value
- `d::Distribution`: A distribution with μ = the weighted mean, and Σ = the weighted variance/covariance matrix
- `blocks_free::Vector{Vector{Int64}}`: A vector of index blocks, where the indices in each block corresponds to the ordering of free parameters only (e.g. all the indices will be ∈ 1:n_free_parameters)
- `blocks_all::Vector{Vector{Int64}}`: A vector of index blocks, where the indices in each block corresponds to the ordering of all parameters (e.g. all the indices will be in ∈ set of all model parameter indices)
- `ϕ_n::Float64`: The current tempering factor
- `ϕ_n1::Float64`: The previous tempering factor

### Keyword Arguments:
- `c::Float64`: The scaling parameter for the proposed covariance matrix
- `α::Float64`: The mixing proportion
- `old_data::Matrix{Float64}`: The matrix of old data to be used in calculating the old_loglh, old_logpost in time tempering

### Outputs:

- `p::Particle`: An updated particle containing updated parameter values, log-likelihood, posterior, and acceptance indicator.

"""
function mutation(m::AbstractModel, data::Matrix{Float64}, p::Particle, d::Distribution,
                  blocks_free::Vector{Vector{Int64}}, blocks_all::Vector{Vector{Int64}},
                  ϕ_n::Float64, ϕ_n1::Float64; c::Float64 = 1., α::Float64 = 1.,
                  old_data::Matrix{Float64} = Matrix{Float64}(size(data, 1), 0))

    n_steps = get_setting(m, :n_mh_steps_smc)

    # draw initial step probability
    # conditions for testing purposes
    step_prob = rand()

    para = p.value
    like = p.loglh
    post_init = p.logpost
    like_prev = p.old_loglh # The likelihood evaluated at the old data (for time tempering)

    # Previous posterior needs to be updated (due to tempering)
    post = post_init + (ϕ_n - ϕ_n1) * like
    accept = false

    for step in 1:n_steps
        for (block_f, block_a) in zip(blocks_free, blocks_all)

            # Index out the parameters corresponding to a given random block
            # And also create a distribution centered at the weighted mean
            # and with Σ corresponding to the same random block
            para_subset = para[block_a]
            d_subset = MvNormal(d.μ[block_f], d.Σ.mat[block_f, block_f])

            para_draw, para_new_density, para_old_density = mvnormal_mixture_draw(para_subset, d_subset;
                                                                                  cc = c, α = α)

            para_new = copy(para)
            para_new[block_a] = para_draw

            like_new = -Inf
            post_new = -Inf
            like_old_data = -Inf
            try
                update!(m, para_new)
                like_new = likelihood(m, data; sampler = true)
                post_new = ϕ_n*like_new + prior(m) - para_new_density
                like_old_data = isempty(old_data) ? 0. : likelihood(m, old_data; sampler = true)
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
    end
    update_mutation!(p, para, like, post, like_prev, accept)
    return p
end
