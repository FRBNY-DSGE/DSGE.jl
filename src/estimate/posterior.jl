"""
`prior(m::AbstractModel{T})`

Calculates log joint prior density of m.parameters.
"""
function prior{T<:AbstractFloat}(m::AbstractModel{T})
    free_params = Base.filter(θ -> !θ.fixed, m.parameters)
    logpdfs = map(logpdf, free_params)
    return sum(logpdfs)
end

"""
```
posterior{T<:AbstractFloat}(m::AbstractModel{T}, data::Matrix{T};
                            mh::Bool = false, catch_errors::Bool = false)
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
function posterior{T<:AbstractFloat}(m::AbstractModel{T},
                                     data::Matrix{T};
                                     mh::Bool = false,
                                     catch_errors::Bool = false)
    catch_errors = catch_errors || mh
    post = likelihood(m, data; mh = mh, catch_errors = catch_errors) + prior(m)
    return post
end

"""
```
posterior!{T<:AbstractFloat}(m::AbstractModel{T}, parameters::Vector{T}, data::Matrix{T};
                             mh::Bool = false, catch_errors::Bool = false)
```

Evaluates the log posterior density at `parameters`.

### Arguments

- `m`: The model object
- `parameters`: New values for the model parameters
- `data`: Matrix of input data for observables

### Optional Arguments

- `mh`: Whether `metropolis_hastings` is the caller. If `mh = true`,
  `catch_errors` is set to `true` (see below)
- `catch_errors`: Whether or not to catch errors of type `GensysError` and
  `ParamBoundsError`
"""
function posterior!{T<:AbstractFloat}(m::AbstractModel{T},
                                      parameters::Vector{T},
                                      data::Matrix{T};
                                      mh::Bool = false,
                                      catch_errors::Bool = false)
    catch_errors = catch_errors || mh
    if mh
        try
            update!(m, parameters)
        catch err
            if isa(err, ParamBoundsError)
                return -Inf
            else
                throw(err)
            end
        end
    else
        update!(m, parameters)
    end
    return posterior(m, data; mh=mh, catch_errors=catch_errors)

end

"""
```
likelihood{T<:AbstractFloat}(m::AbstractModel, data::Matrix{T};
                             mh::Bool = false, catch_errors::Bool = false)
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
function likelihood{T<:AbstractFloat}(m::AbstractModel,
                                      data::Matrix{T};
                                      mh::Bool = false,
                                      catch_errors::Bool = false)
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

    # Return total log-likelihood, excluding the presample
    kal = filter(m, data; catch_errors = catch_errors, allout = false, include_presample = false)

    return kal[:L]
end