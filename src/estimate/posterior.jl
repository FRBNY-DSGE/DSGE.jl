"""
`prior(m::AbstractModel{T})`

Calculates log joint prior density of m.parameters.
"""
function prior(m::AbstractModel{T}) where {T<:AbstractFloat}
    free_params = Base.filter(θ -> !θ.fixed, m.parameters)
    logpdfs = map(logpdf, free_params)
    return sum(logpdfs)
end

"""
```
posterior(m::AbstractModel{T}, data::AbstractArray;
          mh::Bool = false, catch_errors::Bool = false) where {T<:AbstractFloat}
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

- `mh`: Whether `metropolis_hastings` is the caller. If `mh = true`,
  `catch_errors` is set to `true` (see below)
- `catch_errors`: Whether or not to catch errors of type `GensysError` and
  `ParamBoundsError`
"""
function posterior(m::AbstractModel{T},
                   data::AbstractArray;
                   mh::Bool = false,
                   catch_errors::Bool = false) where {T<:AbstractFloat}
    catch_errors = catch_errors || mh
    post = likelihood(m, data; mh = mh, catch_errors = catch_errors) + prior(m)
    return post
end

"""
```
posterior!(m::AbstractModel{T}, parameters::Vector{T}, data::AbstractArray;
           mh::Bool = false, catch_errors::Bool = false) where {T<:AbstractFloat}
```

Evaluates the log posterior density at `parameters`.

### Arguments

- `m`: The model object
- `parameters`: New values for the model parameters
- `data`: AbstractArray of input data for observables

### Optional Arguments

- `mh`: Whether `metropolis_hastings` is the caller. If `mh = true`,
  `catch_errors` is set to `true` (see below)
- `catch_errors`: Whether or not to catch errors of type `GensysError` and
  `ParamBoundsError`
"""
function posterior!(m::AbstractModel{T},
                    parameters::Vector{T},
                    data::AbstractArray;
                    mh::Bool = false,
                    catch_errors::Bool = false) where {T<:AbstractFloat}
    catch_errors = catch_errors || mh
    if mh
        try
            DSGE.update!(m, parameters)
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
    return posterior(m, data; mh=mh, catch_errors=catch_errors)

end

"""
```
likelihood(m::AbstractModel, data::AbstractArray,
           mh::Bool = false, catch_errors::Bool = false) where {T<:AbstractFloat}
```

Evaluate the DSGE likelihood function. Can handle two-part estimation where the observed
sample contains both a normal stretch of time (in which interest rates are positive) and
a stretch of time in which interest rates reach the zero lower bound. If there is a
zero-lower-bound period, then we filter over the 2 periods separately. Otherwise, we
filter over the main sample all at once.

### Arguments

- `m`: The model object
- `data`: matrix of data for observables

### Optional Arguments

- `mh`: Whether `metropolis_hastings` is the caller. If `mh = true`,
  `catch_errors` is set to `true` (see below)
- `catch_errors`: Whether or not to catch errors of type `GensysError`
"""
function likelihood(m::AbstractModel,
                    data::AbstractArray;
                    mh::Bool = false,
                    catch_errors::Bool = false) where {T<:AbstractFloat}
    catch_errors = catch_errors || mh

    # During Metropolis-Hastings, return -∞ if any parameters are not within their bounds
    if mh
        for θ in m.parameters
            (left, right) = θ.valuebounds
            if !θ.fixed && !(left <= θ.value <= right)
                return -Inf
            end
        end
    end

    # Compute state-space system
    system = try
        compute_system(m)
    catch err
        if catch_errors && isa(err, GensysError)
            return -Inf
        else
            rethrow(err)
        end
    end

    # Return total log-likelihood, excluding the presample
    try
        kal = filter(m, data, system; outputs = [:loglh], include_presample = false)
        return kal[:total_loglh]
    catch err
        if catch_errors && isa(err, DomainError)
            @warn "Log of incremental likelihood is negative; returning -Inf"
            return -Inf
        else
            rethrow(err)
        end
    end
end
