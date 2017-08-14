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
                             sampler::Bool = false, catch_errors::Bool = false, φ_smc = 1)
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
-`sampler`: Whether metropolis_hastings or smc is the caller. If `sampler=true`, the log likelihood and the
  transition matrices for the zero-lower-bound period are also returned.
-`catch_errors`: Whether or not to catch errors of type `GensysError` or `ParamBoundsError`
- `φ_smc`: a tempering factor to change the relative weighting of the prior and the
    likelihood when calculating the posterior. It is used primarily in SMC.
"""
function posterior{T<:AbstractFloat}(m::AbstractModel{T},
                                     data::Matrix{T};
                                     sampler::Bool = false,
                                     φ_smc::Float64 = 1.,
                                     catch_errors::Bool = false)
    catch_errors = catch_errors | sampler
    like = likelihood(m, data; sampler=sampler, catch_errors=catch_errors)
    post = φ_smc*like + prior(m)
    if sampler
        return post
    else
        return post
    end
end

"""
```
posterior!{T<:AbstractFloat}(m::AbstractModel{T}, parameters::Vector{T}, data::Matrix{T};
                                  sampler::Bool = false, catch_errors::Bool = false, φ_smc = 1)
```

Evaluates the log posterior density at `parameters`.

### Arguments

- `m`: The model object
- `parameters`: New values for the model parameters
- `data`: Matrix of input data for observables

### Optional Arguments
- `sampler`: Whether metropolis_hastings or smc is the caller. If `sampler=true`, the log likelihood and the
  transition matrices for the zero-lower-bound period are also returned.
- `catch_errors`: Whether or not to catch errors of type `GensysError` or `ParamBoundsError`.
  If `sampler = true`, both should always be caught.
- `φ_smc`: a tempering factor to change the relative weighting of the prior and the
    likelihood when calculating the posterior. It is used primarily in SMC.
"""
function posterior!{T<:AbstractFloat}(m::AbstractModel{T},
                                      parameters::Vector{T},
                                      data::Matrix{T};
                                      sampler::Bool = false,
                                      φ_smc::Float64 = 1.,
                                      catch_errors::Bool = false)
    catch_errors = catch_errors | sampler
    if sampler
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
    return posterior(m, data; sampler=sampler, φ_smc = φ_smc, catch_errors=catch_errors)

end

"""
```
likelihood{T<:AbstractFloat}(m::AbstractModel, data::Matrix{T};
                              sampler::Bool = false, catch_errors::Bool = false)
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
- `sampler`: Whether metropolis_hastings or smc is the caller. If `sampler=true`, the transition matrices for
  the zero-lower-bound period are returned in a dictionary.
- `catch_errors`: If `sampler = true`, `GensysErrors` should always be caught.
"""
function likelihood{T<:AbstractFloat}(m::AbstractModel,
                                      data::Matrix{T};
                                      sampler::Bool = false,
                                      catch_errors::Bool = false)
    catch_errors = catch_errors | sampler

    # During Metropolis-Hastings, return -∞ if any parameters are not within their bounds
    if sampler
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
