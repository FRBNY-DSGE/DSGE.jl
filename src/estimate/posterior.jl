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
- `m`: the model object
- `data`: matrix of data for observables

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
- `m`: The model object
- `parameters`: New values for the model parameters
- `data`: Matrix of input data for observables

### Optional Arguments
- `mh`: Whether metropolis_hastings is the caller. If `mh=true`, the log likelihood and the
  transition matrices for the zero-lower-bound period are also returned.
- `catch_errors`: Whether or not to catch errors of type `GensysError` or `ParamBoundsError`.
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
- `m`: The model object
- `data`: matrix of data for observables

### Optional Arguments
- `mh`: Whether metropolis_hastings is the caller. If `mh=true`, the transition matrices for
  the zero-lower-bound period are returned in a dictionary.
- `catch_errors`: If `mh = true`, `GensysErrors` should always be caught.
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

    # Partition sample into three regimes, and store associated matrices:
    # - R1: presample
    # - R2: normal
    # - R3: zero lower bound and beyond
    R1 = Dict{Symbol, Matrix{T}}()
    R2 = Dict{Symbol, Matrix{T}}()
    R3 = Dict{Symbol, Matrix{T}}()
    regime_mats = [R1, R2, R3]
    regime_likes = zeros(T, 3)

    n_T0         = n_presample_periods(m)
    n_ant        = n_anticipated_shocks(m)
    t_zlb_start  = zlb_start_index(m)
    n_obs_no_ant = n_observables(m) - n_anticipated_shocks(m)
    n_obs        = n_observables(m)
    n_exo        = n_shocks_exogenous(m)

    n_states_no_ant = n_states_augmented(m) - n_anticipated_shocks(m)
    n_states_aug    = n_states_augmented(m)
    n_states        = DSGE.n_states(m)
    regime_states   = [n_states_no_ant, n_states_no_ant, n_states_aug]

    R1[:data] = data[1:n_T0, 1:n_obs_no_ant]
    R2[:data] = data[(n_T0+1):t_zlb_start-1, 1:n_obs_no_ant]
    R3[:data] = data[t_zlb_start:end, :]

    # Step 1: solution to DSGE model - delivers transition equation for the state variables
    # transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
    # If we are in MH, then any errors coming out of gensys should be caught and a -Inf
    # posterior should be returned.
    try
        R3[:TTT], R3[:RRR], R3[:CCC] = solve(m)
    catch err
        if catch_errors && isa(err, GensysError)
            info(err.msg)
            return LIKE_NULL_OUTPUT
        else
            rethrow(err)
        end
    end

    # Get normal, no ZLB matrices
    state_inds = [1:(n_states-n_ant); (n_states+1):n_states_aug]
    shock_inds = 1:(n_exo-n_ant)

    R2[:TTT] = R3[:TTT][state_inds, state_inds]
    R2[:RRR] = R3[:RRR][state_inds, shock_inds]
    R2[:CCC] = R3[:CCC][state_inds, :]

    ## step 2: define the measurement equation: X_t = ZZ*S_t + D + u_t
    ## where u_t = eta_t+MM* eps_t with var(eta_t) = EE
    ## where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

    # Get measurement equation matrices set up for normal and zlb periods
    measurement_R2 = measurement(m, R2[:TTT], R2[:RRR], R2[:CCC]; shocks=false)
    measurement_R3 = measurement(m, R3[:TTT], R3[:RRR], R3[:CCC]; shocks=true)
    for d in (:ZZ, :DD, :QQ, :VVall)
        R2[d] = measurement_R2[d]
        R3[d] = measurement_R3[d]
    end

    # Presample measurement & transition equation matrices are same as normal period
    for d in (:TTT, :RRR, :QQ, :ZZ, :DD, :VVall)
        R1[d] = R2[d]
    end

    ## step 3: compute log-likelihood using Kalman filter
    ##         note that kalman_filter function assumes a transition equation written as:
    ##         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t
    ##         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
    ##         and  VV2 = cov(eps2_t,u_u) = RRR*VV
    ##         define VVall as the joint variance of the two shocks VVall = var([eps2_tu_t])

    # Run Kalman filter on presample
    R1[:A0]         = zeros(T, n_states_no_ant, 1)
    R1[:P0]         = solve_discrete_lyapunov(R1[:TTT], R1[:RRR]*R1[:QQ]*R1[:RRR]')
    out             = kalman_filter(R1[:data]', 1, zeros(T, n_states_no_ant, 1), R1[:TTT], R1[:DD], R1[:ZZ], R1[:VVall], R1[:A0], R1[:P0])
    regime_likes[1] = out[:L]
    R1[:zend]       = out[:zend]
    R1[:Pend]       = out[:Pend]

    # Run Kalman filter on normal period
    zprev           = R1[:zend]
    Pprev           = R1[:Pend]
    out             = kalman_filter(R2[:data]', 1, zeros(regime_states[2], 1), R2[:TTT], R2[:DD], R2[:ZZ], R2[:VVall], zprev, Pprev)
    regime_likes[2] = out[:L]
    R2[:zend]       = out[:zend]
    R2[:Pend]       = out[:Pend]

    # Run Kalman filter on ZLB period
    # This section expands the number of states to accomodate extra states for the
    # anticipated policy shocks. It does so by taking the zend and Pend for the
    # state space without anticipated policy shocks, then shoves in nant
    # zeros in the middle of zend and Pend in the location of
    # the anticipated shock entries.
    before_shocks    = 1:(n_states-n_ant)
    after_shocks_old = (n_states-n_ant+1):(n_states_aug-n_ant)
    after_shocks_new = (n_states+1):n_states_aug

    zprev = [R2[:zend][before_shocks, :];
             zeros(T, n_ant, 1);
             R2[:zend][after_shocks_old, :]]

    Pprev                                     = zeros(T, n_states_aug, n_states_aug)
    Pprev[before_shocks, before_shocks]       = R2[:Pend][before_shocks, before_shocks]
    Pprev[before_shocks, after_shocks_new]    = R2[:Pend][before_shocks, after_shocks_old]
    Pprev[after_shocks_new, before_shocks]    = R2[:Pend][after_shocks_old, before_shocks]
    Pprev[after_shocks_new, after_shocks_new] = R2[:Pend][after_shocks_old, after_shocks_old]

    out             = kalman_filter(R3[:data]', 1, zeros(regime_states[3], 1), R3[:TTT], R3[:DD], R3[:ZZ], R3[:VVall], zprev, Pprev)
    regime_likes[3] = out[:L]
    R3[:zend]       = out[:zend]
    R3[:Pend]       = out[:Pend]

    # Return total log-likelihood, excluding the presample
    like = regime_likes[2] + regime_likes[3]
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
    mats::Dict{Symbol, Matrix{T}}
end
function Posterior{T<:AbstractFloat}(post::T = -Inf,
                                     like::T = -Inf,
                                     mats::Dict{Symbol,Matrix{T}}    = Dict{Symbol,Matrix{T}}())
    return Posterior{T}(post, like, mats)
end
function Base.getindex(P::Posterior, d::Symbol)
    if d in (:post, :like, :mats)
        return getfield(P, d)
    else
        throw(KeyError(d))
    end
end
