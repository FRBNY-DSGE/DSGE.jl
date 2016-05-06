"""
`prior(m::AbstractModel{T})`

Calculates log joint prior density of m.parameters.
"""
function prior{T<:AbstractFloat}(m::AbstractModel{T})
    x = zero(T)
    for θ in m.parameters
        if !θ.fixed
            x += logpdf(θ)
        end
    end
    return x
end

"""
```
posterior{T<:AbstractFloat}(m::AbstractModel{T}, data::Matrix{T};
                             mh::Bool = false, catch_errors::Bool = false)
```

Calculates and returns the log of the posterior distribution for m.parameters:
```
log posterior = log likelihood + log prior
log Pr(Θ|data)  = log Pr(data|Θ)   + log Pr(Θ)
```

### Arguments
-`m`: the model object
-`data`: matrix of data for observables

### Optional Arguments
-`mh`: Whether metropolis_hastings is the caller. If `mh=true`, the log likelihood and the
  transition matrices for the zero-lower-bound period are also returned.
-`catch_errors`: Whether or not to catch errors of type `GensysError` or `ParamBoundsError`
"""
function posterior{T<:AbstractFloat}(m::AbstractModel{T},
                                     data::Matrix{T};
                                     mh::Bool = false,
                                     catch_errors::Bool = false)
    catch_errors = catch_errors | mh
    like, out = likelihood(m, data; mh=mh, catch_errors=catch_errors)
    post = like + prior(m)
    if mh
        return Posterior(post, like, out)
    else
        return Posterior(post, like)
    end
end

"""
```
posterior!{T<:AbstractFloat}(m::AbstractModel{T}, parameters::Vector{T}, data::Matrix{T};
                              mh::Bool = false, catch_errors::Bool = false)
```

Evaluates the log posterior density at `parameters`.

### Arguments
-`m`: The model object
-`parameters`: New values for the model parameters
- `data`: Matrix of input data for observables

### Optional Arguments
-`mh`: Whether metropolis_hastings is the caller. If `mh=true`, the log likelihood and the
  transition matrices for the zero-lower-bound period are also returned.
-`catch_errors`: Whether or not to catch errors of type `GensysError` or `ParamBoundsError`.
  If `mh = true`, both should always be caught.
"""
function posterior!{T<:AbstractFloat}(m::AbstractModel{T},
                                      parameters::Vector{T},
                                      data::Matrix{T};
                                      mh::Bool = false,
                                      catch_errors::Bool = false)
    catch_errors = catch_errors | mh
    if mh
        try
            update!(m, parameters)
        catch err
            return Posterior()
        end
    else
        update!(m, parameters)
    end
    return posterior(m, data; mh=mh, catch_errors=catch_errors)

end

# Empty outputs from `likelihood`.
const LIKE_NULL_DICT   = Dict{Symbol, Matrix{AbstractFloat}}()
const LIKE_NULL_OUTPUT = (-Inf, LIKE_NULL_DICT)

"""
```
likelihood{T<:AbstractFloat}(m::AbstractModel, data::Matrix{T};
                              mh::Bool = false, catch_errors::Bool = false)
```

Evaluate the DSGE likelihood function. Can handle "two part" estimation where the observed
sample contains both a normal stretch of time (in which interest rates are positive) and
a stretch of time in which interest rates reach the zero lower bound. If there is a
zero-lower-bound period, then we filter over the 2 periods separately.  Otherwise, we
filter over the main sample all at once.

### Arguments
-`m`: The model object
-`data`: matrix of data for observables

### Optional Arguments
-`mh`: Whether metropolis_hastings is the caller. If `mh=true`, the transition matrices for
  the zero-lower-bound period are returned in a dictionary.
-`catch_errors`: If `mh = true`, `GensysErrors` should always be caught.
"""
function likelihood{T<:AbstractFloat}(m::AbstractModel,
                                      data::Matrix{T};
                                      mh::Bool = false,
                                      catch_errors::Bool = false)
    catch_errors = catch_errors | mh

    # During Metropolis-Hastings, return -∞ if any parameters are not within their bounds
    if mh
        for θ in m.parameters
            (left, right) = θ.valuebounds
            if !θ.fixed && !(left <= θ.value <= right)
                return LIKE_NULL_OUTPUT
            end
        end
    end

    # Return total log-likelihood, excluding the presample
    R2, R3 = kalman_filter_2part(m, data, allout=false)
    like = R2[:like][1] + R3[:like][1]  # these are matrices with 1 element that we need to extract
    
    if mh
        return like, R3
    else
        return like, LIKE_NULL_DICT
    end
end

# Type returned by posterior
immutable Posterior{T<:AbstractFloat}
    post::T
    like::T
    mats::Dict{Symbol, Array{T}}
end
function Posterior{T<:AbstractFloat}(post::T = -Inf,
                                     like::T = -Inf,
                                     mats::Dict{Symbol,Array{T}}    = Dict{Symbol,Array{T}}())
    return Posterior{T}(post, like, mats)
end
function Base.getindex(P::Posterior, d::Symbol)
    if d in (:post, :like, :mats)
        return getfield(P, d)
    else
        throw(KeyError(d))
    end
end
