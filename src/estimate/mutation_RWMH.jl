"""
```
mutation_RWMH(m, data, p0, l0, post0, i, R, tempering_schedule; rvec=[], rval=[])
```

Execute one proposed move of the Metropolis-Hastings algorithm for a given parameter

### Arguments:
- `p0`: para(j,:)
        initial prior draw.
- `l0`: loglh(j)
        initial log-likelihood at the prior draw.
- `post0`: logpost(j)
           initial log-posterior at the prior draw.
- `tune`: the tune object, which is a struct/type object that contains information about tuning the SMC and the MH algorithms
- `i`: the index of the iteration of the SMC algorithm
- `data`: well-formed data as DataFrame
- `m`: model of type AbstractModel being estimated.

### Keyword Arguments: 
- `rvec`: A matrix of horizontally concatenated random vectors ε as in Θ* = θ + ε for the proposed move in RWMH, for the purposes of testing that mutation_RWMH.
- `rval`: A vector of random values generated in MATLAB for the purposes of testing whether the proposed move exceeds the threshhold and is accepted (or rejected).

### Outputs:

- `ind_para`: para_sim[i,j,:]
              Updated parameter vector
- `ind_loglh`: loglh[j]
               Updated log-likelihood
- `ind_post`: temp_accept[j,1]
              Updated posterior
- `ind_accept`: Indicator for acceptance or rejection (0 reject, 1 accept)

"""
function mutation_RWMH(m::AbstractModel, data::Matrix{Float64}, p0::Array{Float64,1}, l0::Float64, post0::Float64, i::Int64, R::Array{Float64,2}, tempering_schedule::Array{Float64,1}; rvec = [], rval = [])


    n_Φ = get_setting(m, :n_Φ)
    n_part = get_setting(m, :n_particles)
    n_para = n_parameters(m)
    n_blocks = get_setting(m, :n_smc_blocks)

    fixed_para_inds = find([ θ.fixed for θ in m.parameters])
    free_para_inds  = find([!θ.fixed for θ in m.parameters])
    c = get_setting(m, :c)
    accept = get_setting(m, :init_accept)
    target = get_setting(m, :target_accept)

    rvec = isempty(rvec) ? randn(n_para-length(fixed_para_inds),1): rvec
    rval = isempty(rval) ? rand() : rval
    
    cov_mat = chol(R)'
    # SVD is a potential alternative to Cholesky factorization.
    # U, E, V = svd(R)
    # cov_mat = U * diagm(sqrt(E))

    px_draw = p0 + squeeze(c*cov_mat*rvec,2) # reshape to be of the dimension of the original # of paras
    px_draw = reverse(px_draw)
    px = []
    for p in 1:n_para
        if m.parameters[p].fixed == true
            append!(px,[m.parameters[p].value])
        else
            append!(px,[pop!(px_draw)])
        end
    end
    px = convert(Array{Float64},px)
    

    lx = -Inf
    postx = -Inf
    try
        out = posterior!(m,px,data;φ_smc=tempering_schedule[i])
        lx = out[:like]
        postx = out[:post] 
    catch
        # in the event that the proposed move is outside of the bounds or
        # otherwise inappropriate
        lx = -Inf
        postx = -Inf
    end
    
    # Previous posterior needs to be updated (due to tempering)
    post0 = post0+(tempering_schedule[i]-tempering_schedule[i-1])*l0
    

    # Previous prior needs to have the fixed values reattached
    p0_pre = reverse(p0)
    p0 = []
    for p in 1:length(m.parameters)
        if m.parameters[p].fixed == true
            append!(p0,[m.parameters[p].value])
        else
            append!(p0,[pop!(p0_pre)])
        end
    end
    
    # Accept/Reject
    α = exp(postx - post0) # this is RW, so q is canceled out

    if rval .< α # accept
        ind_para   = px 
        ind_loglh  = lx
        ind_post   = postx
        ind_accept = 1
    else
        ind_para   = p0
        ind_loglh  = l0
        ind_post   = post0
        ind_accept = 0
    end

    return ind_para, ind_loglh, ind_post, ind_accept
end
