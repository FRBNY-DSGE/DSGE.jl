"""
```
mutation(m::AbstractModel, data::Matrix{Float64}, p::Particle, R::Array{Float64,2},
         current_φ::Float64, previous_φ::Float64; rvec = [], rval = [], rstep = [])
```

Execute one proposed move of the Metropolis-Hastings algorithm for a given parameter

### Arguments:
- `p_init`: para(j,:)
        initial prior draw.
- `l_init`: loglh(j)
        initial log-likelihood at the prior draw.
- `post_init`: logpost(j)
           initial log-posterior at the prior draw.
- `tune`: the tune object, which is a struct/type object that contains information about tuning the SMC and the MH algorithms
- `i`: the index of the iteration of the SMC algorithm
- `data`: well-formed data as DataFrame
- `m`: model of type AbstractModel being estimated.

### Keyword Arguments:
- `rvec`: A matrix of horizontally concatenated random vectors ε as in Θ* = θ + ε for the proposed move in RWMH, for the purposes of testing that mutation_RWMH.
- `rval`: A vector of random values generated in MATLAB for the purposes of testing whether the proposed move exceeds the threshhold and is accepted (or rejected).

### Outputs:

- `para`: para_sim[i,j,:]
              Updated parameter vector
- `loglh`: loglh[j]
               Updated log-likelihood
- `post`: temp_accept[j,1]
              Updated posterior
- `accept`: Indicator for acceptance or rejection (0 reject, 1 accept)

"""
function mutation(m::AbstractModel, data::Matrix{Float64}, p::Particle, R::Array{Float64,2},
                  current_φ::Float64, previous_φ::Float64; rvec = [], rval = [], rstep = [])

    n_Φ = get_setting(m, :n_Φ)
    n_part = get_setting(m, :n_particles)
    n_para = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_steps = get_setting(m, :n_MH_steps_smc)

    fixed_para_inds = find([ θ.fixed for θ in m.parameters])
    c = get_setting(m, :step_size_smc)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)

    # draw initial step probability
    # conditions for testing purposes
    step_prob = isempty(rval) ? rand() : rval

    U, E, V = svd(R)
    cov_mat = U * diagm(sqrt.(E))

    para = p.value
    like = p.loglh
    post_init = p.logpost

    # Previous posterior needs to be updated (due to tempering)
    post = post_init + (current_φ - previous_φ) * like
    accept = false

    ### BLOCKING ###

    # Generate random ordering
    if m.testing
        para_inds = collect(1:n_para)
    else
        # Returns indices sorted by their random value
        para_inds = [rand() for j in 1:n_para]
        para_inds = sortperm(para_inds)
    end

    # Find the indices for the lower and upper bounds on each block
    b = 1
    blo = Array{Int}(n_blocks)
    bhi = Array{Int}(n_blocks)
    for b in 1:n_blocks
        blo[b] = convert(Int,max(0,ceil((b-1)*n_para/n_blocks)))
        bhi[b] = convert(Int,min(n_para+1,ceil(b*n_para/n_blocks)+1))
    end

    l = 0
    while l < n_steps
        for (j, k) in zip(blo, bhi)
            #Zeroing out relevant rows in cov matrix
            cov_mat_block = copy(cov_mat)
            blo_inds = para_inds[1:j]
            bhi_inds = para_inds[k:end]

            cov_mat_block[:,blo_inds] = 0
            cov_mat_block[blo_inds,:] = 0
            cov_mat_block[:,bhi_inds] = 0
            cov_mat_block[bhi_inds,:] = 0

            para_new =
            if isempty(rstep) && !isempty(rvec)
                para + c*cov_mat_block*rvec
            elseif !isempty(rstep) && isempty(rvec)
                para + c*cov_mat_block*rstep
            else
                rand(DegenerateMvNormal(para, cov_mat_block); cc = c)
            end

            like_new = -Inf
            post_new = -Inf
            try
                update!(m, para_new)
                like_new = likelihood(m, data;sampler = true)
                post_new = current_φ*like_new + prior(m)
            catch err
                if isa(err, GensysError) || isa(err, ParamBoundsError)
                    post_new = like_new = -Inf
                else
                    throw(err)
                end
            end

            # if step is invalid, retry
            if !isfinite(post_new) && isempty(rstep)
                para_new = rand(DegenerateMvNormal(para,cov_mat_block);cc=c)
                step_prob = rand()
                continue
            end

            # Accept/Reject
            α = isempty(rstep) ? exp(post_new - post): 1 # this is RW, so q is canceled out

            if step_prob < α # accept
                para   = para_new
                like   = like_new
                post   = post_new
                accept = true
            end

            # draw again for the next step
            step_prob = rand()
            if m.testing && !isempty(rstep)
                rstep += 1
            end
        end
        l += 1
    end
    update_mutation!(p, para, like, post, accept)
    return p
    # return p,para,like,post,accept
end
