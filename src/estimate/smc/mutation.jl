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
                  block_vars = Vector{Matrix{Float64}}(0, 0)) # RECA: need to add n_mh if want > 1

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
    post = post_init #+ (ϕ_n - ϕ_n1) * like
    accept = 0.0 #false

    for step in 1:n_steps
        eps_step = eps[step, :]
        for (block_f, block_a) in zip(blocks_free, blocks_all)
            step_prob = stepprobs[mm]
            mix_draw  = mixr[mm]
            bvar = block_vars[mm]

            # Index out the parameters corresponding to a given random block
            # And also create a distribution centered at the weighted mean
            # and with Σ corresponding to the same random block
            para_subset = para[block_a]
            d_subset = MvNormal(d.μ[block_f], d.Σ.mat[block_f, block_f])
            eps_block = eps_step[block_a] # RECA

            para_draw, para_new_density, para_old_density = mvnormal_mixture_draw(para_subset,
                                                                                  d_subset;
                                                                                  mixr = mix_draw, # RECA
                                                                                  eps = eps_block, # RECA
                                                                                  c=c, α=α, mu=mu[block_a],
                                                                                  bvar=bvar, pnum=pnum)
            para_new = deepcopy(para)
            para_new[block_a] = para_draw

            ##### BEGINING OF NONSENSE #######
            lik0 = like
            pr0  = post

            q0 = α * exp(logpdf(MvNormal(para_draw,   c^2 * d_subset.Σ.mat), para_subset))
            q1 = α * exp(logpdf(MvNormal(para_subset, c^2 * d_subset.Σ.mat), para_draw))

            if (pnum ==  12001) @show "1. ", α, q0, q1 end

            ind_pdf = 1.
            for i=1:length(block_a)
                sigi = sqrt(d_subset.Σ.mat[i,i])

                if (pnum ==  12001) @show "2.5 ", i, para_subset[i], para_draw[i] end

                zstat = (para_subset[i] - para_draw[i]) / sigi
                ind_pdf = ind_pdf / (sigi * sqrt(2. * 3.1415)) * exp(-0.5 * zstat ^ 2)

                if (pnum ==  12001) @show "2. ", i, sigi, zstat, ind_pdf end
            end

            q0 = q0 + (1. - α) / 2. * ind_pdf
            q1 = q1 + (1. - α) / 2. * ind_pdf

            if (pnum ==  12001) @show "3. ", q0, q1 end

            q0 = q0 + (1. - α) / 2. * exp(logpdf(MvNormal(mu[block_a], c^2 * d_subset.Σ.mat), para_subset))
            q1 = q1 + (1. - α) / 2. * exp(logpdf(MvNormal(mu[block_a], c^2 * d_subset.Σ.mat), para_draw))

            if (pnum ==  12001) @show "4. ", pnum, q0, q1 end

            q0 = log(q0)
            q1 = log(q1)

            if (pnum ==  12001) @show "5. ", q0, q1 end

            lik1 = -Inf
            pr1  = -Inf

            n_para = length(para)
            if (pnum == 12001)
                @show "TEST", "BLOCKING BELOW"
                p1 = readdlm("/home/rcerxs30/SLICOT-2018-12-19/dsge-smc/fortran/smc-sw-new-mix-npart-12000-nintmh-1-nphi-500-prior-b3-trial1-phibend-jan28-mixron-blocking-FIX/002step_p1.txt")
                fort_para = vec(p1[(mm-1)*n_steps*n_para + (step-1)*n_para + 1:(mm-1)*n_steps*n_para + step*n_para])
                update!(m, fort_para)
                @show likelihood(m, data)
            end
            try
                update!(m, para_new)
                pr1 = prior(m)
                lik1 = likelihood(m, data; sampler=false, use_chand_recursion=false, verbose=verbose)
                if (pnum == 12001) @show lik1 end
                if (pnum == 12001) @show isapprox(fort_para, para_new, atol=1e-4) end
            catch
                pr1 = lik1 = -Inf
            end
            if (q0 == Inf && q1 == Inf)
                q0 = 0.
            end
            if (pnum ==  12001) @show "6. ", pr1, lik1 end
            if (pnum ==  12001) @show ϕ_n, lik1, lik0, pr1, pr0, q0, q1 end
            if (pnum ==  12001) @show ϕ_n * (lik1 - lik0) + pr1 - pr0 + q0 - q1 end
            if (pnum ==  12001) @show ϕ_n * (lik1 - lik0) end
            if (pnum ==  12001) @show pr1, pr0, pr1 - pr0 end
            if (pnum ==  12001) @show q0, q1, q0 - q1 end

            my_alp = exp(ϕ_n * (lik1 - lik0) + pr1 - pr0 + q0 - q1)

            if my_alp > 1.0
                my_alp = 1.0
            elseif isnan(my_alp)
                my_alp = 0.0
            end
            #@show "7. ", my_alp, alp[1]
            #@show pnum, mm, my_alp, alp[mm], pr0
###################### old code goes here: old_mutation_body.jl
            if !(abs(alp[mm] - my_alp) < 1e-3)
                @show pnum, mm, alp[mm], my_alp
                @show q1, q0, lik1, lik0, pr1, pr0, step_lik[mm]
            end

            if abs(alp[mm] - my_alp) / abs(alp[mm]) > 1e-3 && abs(alp[mm] - my_alp) > 1e-3
                @show pnum, mm, alp[mm], my_alp
            end
            if alp[mm] > 0.1
               #@assert abs(alp[mm] - my_alp) / abs(alp[mm]) < 5e-2
                if !(abs(alp[mm] - my_alp) / abs(alp[mm]) < 5e-2)
                    print_with_color(:red, "Assert Fails abs(alp[mm] - my_alp) / abs(alp[mm]) < 5e-2.\n")
                    #@show pnum
                    #error()
                end
            else
               #@assert abs(alp[mm] - my_alp) < 1.5e-2
                if !(abs(alp[mm] - my_alp) < 1.5e-2)
                    print_with_color(:red, "Assert Fails abs(alp[mm] - my_alp) < 1.5e-2.\n")
                end
            end

            if step_prob < alp[mm]#my_alp# JUST ALP running in top pane
                para = para_new
                like = lik1 # BOTH ALP and resetting lik in middle pane
                post = pr1
                accept += length(block_a)
            end

            #@show "8. ", para, like, post

            # draw again for the next step
            #step_prob = rand()

            if !isapprox(step_pr[(mm-1)*n_steps*n_para + (step-1)*n_para + 1:(mm-1)*n_steps*n_para + step*n_para], para, atol=1e-3)
                @show step_pr[(mm-1)*n_steps*n_para + (step-1)*n_para + 1:(mm-1)*n_steps*n_para + step*n_para], para
                @show abs.(step_pr[(mm-1)*n_steps*n_para + (step-1)*n_para + 1:(mm-1)*n_steps*n_para + step*n_para] - para)
                print_with_color(:red, "Assert Fails on ("*string(pnum)*", "*string(mm)*"): step_pr ≈ para.\n")
                #@assert step_pr ≈ para
            end

            #if !(abs(like - step_lik[mm]) / abs(like) < 1e-3)
            #    @show pnum, q1, q0, lik1, lik0, pr1, pr0, alp[mm], my_alp, step_prob
            #    @show pnum, abs(like - step_lik[mm]) / abs(like), like, step_lik[mm]
            #end
            if (like <= 1e-3 && step_lik[mm] <= 1e-3)
                if !(abs(like - step_lik[mm]) / abs(like) < 1e-3)
                    @show "BAD NEWS", pnum, abs(like - step_lik[mm]) / abs(like), like, step_lik[mm], alp[mm], my_alp, step_prob
                else
                    if !(abs(like - step_lik[mm]) / abs(like) < 1e-3)
                        print_with_color(:red, "Assert Fails abs(like - step_lik[mm]) / abs(like) < 1e-3.\n")
                    end
                end
            else
                if !(abs(like - step_lik[mm]) / abs(like) < 1e-3)
                    print_with_color(:red, "Assert Fails abs(like - step_lik[mm]) / abs(like) < 1e-3.\n")
                end
            end
            #like = step_lik[mm]
            mm += 1
        end
    end
    if (pnum % 1000 == 0) @show pnum end
    update_mutation_FORTRAN!(p, para, like, post, like_prev, accept / 36.0) # RECA: hardcoded
    return p
end
