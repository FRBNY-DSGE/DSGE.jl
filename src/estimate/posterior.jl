"""
`prior(m::AbstractDSGEModel{T})`

Calculates log joint prior density of m.parameters.
"""
function prior(m::AbstractDSGEModel{T}) where {T<:AbstractFloat}
    return ModelConstructors.prior(m.parameters)
end

"""
```
posterior(m::Union{AbstractDSGEModel{T},AbstractVARModel{T}}, data::Matrix{T};
          sampler::Bool = false, catch_errors::Bool = false,
          φ_smc::Float64 = 1) where {T<:AbstractFloat}
```

Calculates and returns the log of the posterior distribution for `m.parameters`:

```
log posterior  = log likelihood + log prior + const
log Pr(Θ|data) = log Pr(data|Θ) + log Pr(Θ) + const
```

### Arguments

- `m`: the model object
- `data`: matrix of data for observables

### Optional Arguments
-`sampler`: Whether metropolis_hastings or smc is the caller. If `sampler=true`,
    the log likelihood and the transition matrices for the zero-lower-bound
    period are also returned.
-`catch_errors`: Whether to catch errors of type `GensysError` or `ParamBoundsError`
- `φ_smc`: a tempering factor to change the relative weighting of the prior and
     the likelihood when calculating the posterior. It is used primarily in SMC.
"""
function posterior(m::Union{AbstractDSGEModel{T},AbstractVARModel{T}}, data::AbstractArray;
                   sampler::Bool = false, ϕ_smc::Float64 = 1.,
                   catch_errors::Bool = false) where {T<:AbstractFloat}
    catch_errors = catch_errors | sampler
    like = likelihood(m, data; sampler = sampler,
                      catch_errors = catch_errors)
    post = ϕ_smc*like + prior(m)
    return post
end

"""
```
posterior!(m::Union{AbstractDSGEModel{T},AbstractVARModel{T}},
           parameters::Vector{T}, data::Matrix{T};
           sampler::Bool = false, catch_errors::Bool = false,
           ϕ_smc::Float64 = 1., toggle::Bool = true) where {T<:AbstractFloat}
```

Evaluates the log posterior density at `parameters`.

### Arguments

- `m`: The model object
- `parameters`: New values for the model parameters
- `data`: Matrix of input data for observables

### Optional Arguments
- `sampler`: Whether metropolis_hastings or smc is the caller. If `sampler=true`,
     the log likelihood and the transition matrices for the zero-lower-bound
     period are also returned.
- `catch_errors`: Whether to catch errors of type `GensysError` or `ParamBoundsError`
     If `sampler = true`, both should always be caught.
- `ϕ_smc`: a tempering factor to change the relative weighting of the prior and
     the likelihood when calculating the posterior. It is used primarily in SMC.
- `toggle`: if true, we call `ModelConstructors.toggle_regime!(values)` before
    updating any values to ensure the `value` field of the parameters in `values`
    correspond to regime 1 values.
"""
function posterior!(m::Union{AbstractDSGEModel{T},AbstractVARModel{T}},
                    parameters::Vector{T}, data::AbstractArray;
                    sampler::Bool = false, ϕ_smc::Float64 = 1.,
                    catch_errors::Bool = false, toggle::Bool = true) where {T<:AbstractFloat}
    catch_errors = catch_errors | sampler
    regime_switching = haskey(get_settings(m), :regime_switching) && get_setting(m, :regime_switching)
    if sampler
        try
            DSGE.update!(m, parameters; regime_switching = regime_switching, toggle = toggle)
        catch err
            if isa(err, ParamBoundsError)
                return -Inf
            else
                throw(err)
            end
        end
    else
        DSGE.update!(m, parameters)
    end
    return posterior(m, data; sampler = sampler, ϕ_smc = ϕ_smc, catch_errors = catch_errors)
end

"""
```
likelihood(m::AbstractDSGEModel, data::Matrix{T};
           sampler::Bool = false, catch_errors::Bool = false) where {T<:AbstractFloat}
```

Evaluate the DSGE log-likelihood function. Can handle two-part estimation where the observed
sample contains both a normal stretch of time (in which interest rates are positive) and
a stretch of time in which interest rates reach the zero lower bound. If there is a
zero-lower-bound period, then we filter over the 2 periods separately. Otherwise, we
filter over the main sample all at once.

### Arguments

- `m`: The model object
- `data`: matrix of data for observables

### Optional Arguments
- `sampler`: Whether metropolis_hastings or smc is the caller. If `sampler=true`, the
    transition matrices for the zero-lower-bound period are returned in a dictionary.
- `catch_errors`: If `sampler = true`, `GensysErrors` should always be caught.
"""
function likelihood(m::AbstractDSGEModel, data::AbstractMatrix;
                    sampler::Bool = false,
                    catch_errors::Bool = false,
                    use_chand_recursion::Bool = false,
                    tol::Float64 = 0.0,
                    verbose::Symbol = :high) where {T<:AbstractFloat}

    catch_errors = catch_errors | sampler
    use_penalty  = try get_setting(m, :use_likelihood_penalty) catch; false end
    auto_reject  = try get_setting(m, :auto_reject) catch; false end

    if auto_reject
        m <= Setting(:auto_reject, false)
        return -Inf
    end

    # During Metropolis-Hastings, return -∞ if any parameters are not within their bounds
    if sampler
        for θ in m.parameters
            (left, right) = θ.valuebounds
            if !θ.fixed && !(left <= θ.value <= right)
                return -Inf
            end
        end
    end

    # Likelihood penalties
    ψ_l, ψ_p, penalty = 1.0, 1.0, 0.0
    if use_penalty
        ψ_l         = get_setting(m, :ψ_likelihood)
        ψ_p         = get_setting(m, :ψ_penalty)
        target_vars = get_setting(m, :target_vars)
        target_σt   = get_setting(m, :target_σt)
        targets     = get_setting(m, :targets)

        for (var, target, σt) in zip(target_vars, targets, target_σt)
            try
                penalty +=  -0.5 * (log(target) - log(m[var].value))^2 / σt^2
            catch err
                println(err)
                return -Inf
            end
        end
        if ψ_l == 0.0
            return ψ_p * penalty
        end
    end

    # Compute state-space system
    system = try
        compute_system(m; tvis = haskey(get_settings(m), :tvis_information_set), verbose = verbose)
    catch err
        if catch_errors && (isa(err, GensysError) || isa(err, KleinError))
            return -Inf
        else
            rethrow(err)
        end
    end

    # Return total log-likelihood, excluding the presample
    try
        if isa(m, PoolModel)
            return ψ_l * sum(filter_likelihood(m, data; tol = tol,
                                               tuning = get_setting(m, :tuning))) + ψ_p * penalty
        elseif use_chand_recursion==false
            return ψ_l * sum(filter_likelihood(m, data, system;
                                               include_presample = false, tol = tol)) +
                                                   ψ_p * penalty
        else
            # TODO: extend chand_recursion for a regime-switching system
            return ψ_l * chand_recursion(data, system[:TTT], system[:RRR], system[:CCC],
                                   system[:QQ], system[:ZZ], system[:DD], system[:EE];
                                   allout = true, Nt0 = n_presample_periods(m),
                                   tol = tol)[1] + ψ_p * penalty
        end
    catch err
        if catch_errors && isa(err, DomainError)
            @warn "Log of incremental likelihood is negative; returning -Inf"
            return -Inf
        else
            rethrow(err)
        end
    end
end

"""
```
likelihood(m::AbstractVARModel, data::Matrix{T};
           sampler::Bool = false, catch_errors::Bool = false,
           verbose::Symbol = :high) where {T<:AbstractFloat}
```

Evaluate a VAR likelihood function.

### Arguments
- `m`: The model object
- `data`: matrix of data for observables

### Optional Arguments
- `sampler`: Whether metropolis_hastings or smc is the caller. If `sampler=true`, the
    transition matrices for the zero-lower-bound period are returned in a dictionary.
- `catch_errors`: If `sampler = true`, `GensysErrors` should always be caught.
"""
function likelihood(m::AbstractVARModel, data::AbstractMatrix;
                    sampler::Bool = false,
                    catch_errors::Bool = false,
                    verbose::Symbol = :high) where {T<:AbstractFloat}

    catch_errors = catch_errors | sampler
    use_penalty  = try get_setting(m, :use_likelihood_penalty) catch; false end
    auto_reject  = try get_setting(m, :auto_reject) catch; false end

    if auto_reject
        m <= Setting(:auto_reject, false)
        return -Inf
    end

    # During Metropolis-Hastings, return -∞ if any parameters are not within their bounds
    if sampler
        if isa(m, AbstractDSGEVARModel)
            for θ in get_parameters(m)
                (left, right) = θ.valuebounds
                if !θ.fixed && !(left <= θ.value <= right)
                    return -Inf
                end
            end
        end
    end

    # Likelihood penalties
    ψ_l, ψ_p, penalty = 1.0, 1.0, 0.0
    if use_penalty
        ψ_l         = get_setting(m, :ψ_likelihood)
        ψ_p         = get_setting(m, :ψ_penalty)
        target_vars = get_setting(m, :target_vars)
        target_σt   = get_setting(m, :target_σt)
        targets     = get_setting(m, :targets)

        for (var, target, σt) in zip(target_vars, targets, target_σt)
            try
                penalty +=  -0.5 * (log(target) - log(m[var].value))^2 / σt^2
            catch err
                println(err)
                return -Inf
            end
        end
        if ψ_l == 0.0
            return ψ_p * penalty
        end
    end

    # Return total log-likelihood (presample for VAR is excluded)
    try
        if isa(m, AbstractDSGEVECMModel)
            return ψ_l * dsgevecm_likelihood(m, data) + ψ_p * penalty
        elseif isa(m, AbstractDSGEVARModel)
            return ψ_l * dsgevar_likelihood(m, data) + ψ_p * penalty
        end
    catch err
        if catch_errors && (isa(err, GensysError) || isa(err, KleinError))
            return -Inf
        elseif catch_errors && isa(err, DomainError)
            @warn "Log of incremental likelihood is negative; returning -Inf"
            return -Inf
        else
            rethrow(err)
        end
    end
end
