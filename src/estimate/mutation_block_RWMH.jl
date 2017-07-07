"""
```
mutation_block_RWMH(m, data, p0, l0, post0, i, R, tempering_schedule; rvec=[], rval=[])
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
function mutation_block_RWMH(m::AbstractModel, data::Matrix{Float64}, para_init::Array{Float64,1}, like_init::Float64, post_init::Float64, i::Int64, R::Array{Float64,2}, tempering_schedule::Array{Float64,1}; rvec = [], rval = [])


    n_Φ = get_setting(m, :n_Φ)
    n_part = get_setting(m, :n_particles)
    n_para = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_steps = get_setting(m, :n_MH_steps_smc)

    fixed_para_inds = find([ θ.fixed for θ in m.parameters])
    free_para_inds  = find([!θ.fixed for θ in m.parameters])
    c = get_setting(m, :step_size_smc)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)

    # draw initial step and step probability
    step = isempty(rvec) ? randn(n_para-length(fixed_para_inds),1): rvec
    step_prob = isempty(rval) ? rand() : rval

    cov_mat = chol(R)'

    para = para_init
    like = like_init
    # Previous posterior needs to be updated (due to tempering)
    post = post_init + (tempering_schedule[i] - tempering_schedule[i-1]) * like
    accept = 0
    
### BLOCKING ###

    # Generate random index for each parameter, as means of sorting
    para_idx = zeros(n_para)
    for i=1:n_para
        para_idx[i] = rand()
    end
    # Returns indices sorted by their random value
    para_idx = sortperm(para_idx)
    @show para_idx
    
### BLOCKING ###

    j = 0
    b = 1
    while j < n_steps
                
        while b <= n_blocks # corresponds to jth block
            if b > 1
                cov_mat_temp = isolate_block(cov_mat, b, n_blocks, para_idx)
                @show cov_mat_temp
                para_new = para + c * cov_mat_temp * step
            else
                para_new = para + c * cov_mat * step
            end
            @show para_new
            @show n_blocks
            post_new = posterior!(m, augment_draw(para_new, m), data; 
                                  φ_smc = tempering_schedule[i],
                                  sampler = true)
            like_new = post_new - prior(m)

            # if step is invalid, retry
            if !isfinite(post_new)
                #j -= 1
                b -= 1    
            end

            # Accept/Reject
            α = exp(post_new - post) # this is RW, so q is canceled out

            if step_prob < α # accept
                para   = para_new
                like   = like_new
                post   = post_new
                accept = 1
            end

            # draw again for the next step
            step = randn(n_para-length(fixed_para_inds),1)
            step_prob = rand()
            # j += 1
            b += 1

        end
    end
    return augment_draw(para, m), like, post, accept
end

# Zeros-out the lines of the matrices not in block
function isolate_block(cov_mat::AbstractArray, b, n_blocks, para_idx::AbstractArray)
    cov_mat_block = cov_mat
    n_para = length(para_idx)
    sup = min(b*n_blocks,n_para)
    inf = max(1,(b-1)*n_blocks)
    for i = 1:inf
        idx = para_idx[i]
        cov_mat_block[:,idx] = zeros(n_para)
        cov_mat_block[idx,:] = zeros(n_para)
    end
    for i = sup:n_para
        idx = para_idx[i]
        cov_mat_block[:,idx] = zeros(n_para)
        cov_mat_block[idx,:] = zeros(n_para)
    end
    return cov_mat_block
end


function augment_draw(draw, m)
    new_draw = map(x -> x.value, m.parameters)
    j = 1
    for (i,θ) in enumerate(m.parameters)
        if !θ.fixed
            new_draw[i] = draw[j]
            j += 1
        end
    end
    return new_draw
end
