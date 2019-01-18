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
                  old_data::Matrix{Float64} = Matrix{Float64}(size(data, 1), 0),
                  use_chand_recursion::Bool = false,
                  verbose::Symbol = :low,
                  pnum::Int64 = -1,
                  mixr = Vector{Float64}(size(blocks_all,1), 0), # RECA: need to add n_mh if want > 1
                  eps  = Matrix{Float64}(size(1, 41), 0), # RECA : size(n_steps, n_para)
                  stepprobs = Vector{Float64}(size(blocks_all,1), 0),
                  step_pr = Matrix{Float64}(size(1, 1), 0),
                  step_lik = Matrix{Float64}(size(1, 1), 0),
                  step_p1 = Matrix{Float64}(size(1, 1), 0),
                  alp = Matrix{Float64}(size(1, 1), 0),
                  mu = Vector{Float64}(size(41,1), 0), # RECA: need to add n_mh if want > 1
                  bvar = Matrix{Float64}(0, 0)) # RECA: need to add n_mh if want > 1

    n_steps = get_setting(m, :n_mh_steps_smc)

    # draw initial step probability
    # conditions for testing purposes
    #step_prob = rand()
    mm = 1 # RECA

    para = p.value
    like = p.loglh
    post_init = p.logpost
    like_prev = p.old_loglh # The likelihood evaluated at the old data (for time tempering)

    # Previous posterior needs to be updated (due to tempering)
    post = post_init + (ϕ_n - ϕ_n1) * like
    accept = 0.0 #false

    #println("***************************************************** "*string(pnum))

    for step in 1:n_steps
        eps_step = eps[step, :]
        for (block_f, block_a) in zip(blocks_free, blocks_all)
            step_prob = stepprobs[mm]
            mix_draw  = mixr[mm]

            # Index out the parameters corresponding to a given random block
            # And also create a distribution centered at the weighted mean
            # and with Σ corresponding to the same random block
            para_subset = para[block_a]
            d_subset = MvNormal(d.μ[block_f], d.Σ.mat[block_f, block_f])

            # RECA
            eps_block = eps_step[block_a]
            para_draw, para_new_density, para_old_density = mvnormal_mixture_draw(para_subset,
                                                                                  d_subset;
                                                                                  mixr = mix_draw, # RECA
                                                                                  eps = eps_block, # RECA
                                                                                  cc=c, α=α, mu=mu[block_a],
                                                                                  bvar = bvar)
            para_new = deepcopy(para)
            para_new[block_a] = para_draw
            @assert para_new ≈ step_p1

            ##### BEGINING OF NONSENSE #######
            p0 = para
            lik0 = like
            pr0 = post_init
            function make_sym(mat::Matrix{Float64})
                return (mat + mat') / 2
            end
            q0 = α * exp(logpdf(MvNormal(para_draw, c^2*d_subset.Σ.mat), para_subset)) #make_sym(bvar'*bvar)), para_subset))
            q1 = α * exp(logpdf(MvNormal(para_subset, c^2*d_subset.Σ.mat), para_draw))  #make_sym(bvar'*bvar)), para_draw))

            #@show "1. ", α, q0, q1

            ind_pdf = 1.
            for i=1:length(block_a)
                sigi = sqrt(d_subset.Σ.mat[i,i])
                #@show "2.5 ", i, para_subset[i], para_draw[i]
                zstat = (para_subset[i] - para_draw[i]) / sigi
                ind_pdf = ind_pdf / (sigi * sqrt(2. * 3.1415)) * exp(-0.5 * zstat ^ 2)
                #@show "2. ", i, sigi, zstat, ind_pdf
            end

            q0 = q0 + (1. - α)/2. * ind_pdf
            q1 = q1 + (1. - α)/2. * ind_pdf

            #@show "3. ", q0, q1

            q0 = q0 + (1. - α)/2. * exp(logpdf(MvNormal(mu[block_a], c^2 * d_subset.Σ.mat), para_subset))  #make_sym(bvar'*bvar)), para_subset))
            q1 = q1 + (1. - α)/2. * exp(logpdf(MvNormal(mu[block_a], c^2 * d_subset.Σ.mat), para_draw)) #make_sym(bvar'*bvar)), para_draw))

            #@show "4. ", q0, q1

            q0 = log(q0)
            q1 = log(q1)

            #@show "5. ", q0, q1

            lik1 = -Inf
            pr1  = -Inf

            try
                update!(m, para_new)
                pr1 = prior(m)
                lik1 = likelihood(m, data; sampler=true, use_chand_recursion=use_chand_recursion, verbose=verbose)
            catch
                #@show "SETTING TO -INF"
                pr1 = lik1 = -Inf
            end
            if (q0 == Inf && q1 == Inf)
                q0 = 0.
                #q1 = 0.
            end
            #@show "6. ", pr1, lik1

            #@show ϕ_n, lik1, lik0, pr1, pr0, q0, q1
            #@show ϕ_n * (lik1 - lik0) + pr1 - pr0 + q0 - q1
            #@show ϕ_n * (lik1 - lik0)
            #@show pr1, pr0, pr1 - pr0
            #@show q0, q1, q0 - q1
            my_alp = exp(ϕ_n * (lik1 - lik0) + pr1 - pr0 + q0 - q1)

            if my_alp > 1.0
                my_alp = 1.0
            end
            #@show "7. ", my_alp, alp[1]

            ##### END OF NONSENSE, REAL CODE BELOW #######
#=
            like_new = -Inf
            post_new = -Inf
            like_old_data = -Inf

            try
                update!(m, para_new)
                like_new = likelihood(m, data; sampler = true,
                                      use_chand_recursion = use_chand_recursion,
                                      verbose = verbose)
                if like_new == -Inf
                    post_new = like_old_data = -Inf
                end
                post_new = ϕ_n*like_new + prior(m) - para_new_density
                like_old_data = isempty(old_data) ? 0. : likelihood(m, old_data; sampler = true,
                                                         use_chand_recursion = use_chand_recursion,
                                                         verbose = verbose)
            catch err
                if isa(err, ParamBoundsError)
                    println("ParamBoundsError!")
                    post_new = like_new = like_old_data = -Inf
                else
                    throw(err)
                end
            end

            # Accept Reject
            post_old = post - para_old_density
            η        = exp(post_new - post_old)

            if η > 1.
                η = 1.0
            end
            #@show  η, alp[1]
            #@show step_pr, step_lik
            #@show post_new, post_old, post, para_old_density

            if step_prob < η # accept
                para      = para_new
                like      = like_new
                post      = post_new + para_new_density # Have to add it back so as to not accumulate
                like_prev = like_old_data               # para_new_density throughout the iterations
                accept   += length(block_a) #true
            end

            @show η, alp[1], my_alp
            @assert abs(η - alp[1]) < 1e-4

            @show like_new, lik1
            @show post_new, para_new_density, pr1

            if !(like_new == -Inf && lik1 == -Inf) @assert abs(like_new - lik1) < 1e-4 end
            if !(post_new == -Inf && pr1 == -Inf) @assert abs(post_new + para_new_density - pr1) < 10 end
=#
################################# REAL CODE ABOVE #########
            @show pnum
            @assert abs(alp[1] - my_alp) < 1e-4

            if step_prob < my_alp
                para = para_new
                like = lik1
                post = pr1
                accept += length(block_a)
            end

            #@show "8. ", para, like, post

            mm += 1
            # draw again for the next step
            @assert step_pr ≈ para
            @assert like ≈ step_lik[1]


            #step_prob = rand()
        end
    end
    update_mutation_FORTRAN!(p, para, like, post, like_prev, accept / 41.0) # RECA: hardcoded
    return p
end
