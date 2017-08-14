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
function mutation_block_RWMH(m::AbstractModel, data::Matrix{Float64}, para_init::Array{Float64,1}, like_init::Float64, post_init::Float64, i::Int64, R::Array{Float64,2}, tempering_schedule::Array{Float64,1}; rvec = [], rval = [], rstep = [])

    n_Φ = get_setting(m, :n_Φ)
    n_part = get_setting(m, :n_particles)
    n_para = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)
    n_steps = get_setting(m, :n_MH_steps_smc)

    fixed_para_inds = find([ θ.fixed for θ in m.parameters])
    c = get_setting(m, :step_size_smc)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)

    # draw initial step and step probability
    # conditions for testing purposes
    # step = isempty(rstep) ? (isempty(rvec) ? randn(n_para,1): rvec) : rstep
    step_prob = isempty(rval) ? rand() : rval

    cov_mat = R
    para = para_init
    like = like_init

    # Previous posterior needs to be updated (due to tempering)
    post = post_init + (tempering_schedule[i] - tempering_schedule[i-1]) * like
    accept = 0
    
### BLOCKING ###

    # Generate random index for each parameter, as means of sorting
    if m.testing
        para_inds = collect(1:n_para) 
    else
        para_inds = zeros(n_para)
        for j=1:n_para
            para_inds[j] = rand()
        end
        # Returns indices sorted by their random value
        para_inds = sortperm(para_inds)
    end

    # Find the indices for the lower and upper bounds on each block
    b = 1
    blo = Array{Int64}(n_blocks)
    bhi = Array{Int64}(n_blocks)
    for b in 1:n_blocks
        blo[b] = convert(Int64,max(0,ceil((b-1)*n_para/n_blocks)))
        bhi[b] = convert(Int64,min(n_para+1,ceil(b*n_para/n_blocks)+1))
    end

    for (j,k) in zip(blo,bhi)
        #Zeroing out relevant rows in cov matrix
        cov_mat_block = copy(cov_mat)
        blo_inds = para_inds[1:j]
        bhi_inds = para_inds[k:end]
        cov_mat_block[:,blo_inds] = zeros(size(cov_mat_block[:,blo_inds]))
        cov_mat_block[blo_inds,:] = zeros(size(cov_mat_block[blo_inds,:]))
        cov_mat_block[:,bhi_inds] = zeros(size(cov_mat_block[:,bhi_inds]))
        cov_mat_block[bhi_inds,:] = zeros(size(cov_mat_block[bhi_inds,:]))

        l = 0
        while l < n_steps
            #para_new = para+c*cov_mat_block*step
            para_new = rand(DegenerateMvNormal(para,cov_mat_block);cc=c)
            post_new = posterior!(m, vec(para_new), data; 
                                  φ_smc = tempering_schedule[i],
                                  sampler = true)
            like_new = post_new - prior(m)

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
                accept = 1
            end

            # draw again for the next step
            #step = isempty(rstep) ? randn(n_para,1): step.+1
            step_prob = rand()
            l += 1
        end
    end
    return para, like, post, accept
end

# Function used at the start of program to zero-out fixed rows and columns
function augment_cov_mat(cov_mat,list_inds)
    for id in list_inds
        cov_mat[:,id] = zeros(size(cov_mat)[1])
        cov_mat[id,:] = zeros(size(cov_mat)[1])
    end
    return cov_mat
end

# Not currently being called
function augment_draw(draw, m)
    new_draw = map(x -> x.value, m.parameters)
    #j = 1
    #for (i,θ) in enumerate(m.parameters)
    #    if !θ.fixed
    #        new_draw[i] = draw[j]
    #        j += 1
    #    end
    #end
    return new_draw
end
