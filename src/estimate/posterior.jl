#=
doc"""
Calculates (log of) the joint density of the model parameters.
"""
=#
function prior{T<:AbstractFloat}(model::AbstractDSGEModel{T})
    x = zero(T)
    for θ in model.parameters
        if !θ.fixed
            x += logpdf(θ)
        end
    end
    return x
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

#=
doc"""
posterior{T<:AbstractFloat}(model::AbstractDSGEModel, YY::Matrix{T}; mh::Bool = false, catch_errors::Bool = false)

### Parameters
-`model`: the model object
-`YY`: matrix of data for observables

### Optional Arguments
-`mh`: Whether metropolis_hastings is the caller. If `mh=true`, the log likelihood and the transition matrices for the zero-lower-bound period are also returned.
-`catch_errors`: Whether or not to catch errors of type `GensysError` or `ParamBoundsError`

### Description
Calculates and returns the log of the posterior distribution for the model parameters:
  log posterior = log likelihood + log prior
  log Pr(Θ|YY)  = log Pr(YY|Θ)   + log Pr(Θ)    # where Θ is `m.parameters`
"""
=#
function posterior{T<:AbstractFloat}(model::AbstractDSGEModel,
                                     YY::Matrix{T};
                                     mh::Bool = false,
                                     catch_errors::Bool = false)
    catch_errors = catch_errors | mh
    like, out = likelihood(model, YY; mh=mh, catch_errors=catch_errors)
    post = like + prior(model)
    if mh
        return Posterior(post, like, out)
    else
        return Posterior(post, like)
    end
end

#=
doc"""
posterior!{T<:AbstractFloat}(model::AbstractDSGEModel, parameters::Vector{T}, YY::Matrix{T}; [mh::Bool = false], [catch_errors::Bool = false])

### Parameters
-`model`: The model object
-`parameters`: New values for the model parameters
- `YY`: Matrix of input data for observables

### Optional Arguments
-`mh`: Whether metropolis_hastings is the caller. If `mh=true`, the log likelihood and the transition matrices for the zero-lower-bound period are also returned.
-`catch_errors`: Whether or not to catch errors of type `GensysError` or `ParamBoundsError`. If `mh = true`, both should always be caught.

### Description
Evaluates the log posterior distribution at `parameters`
"""
=#
function posterior!{T<:AbstractFloat}(model::AbstractDSGEModel{T},
                                      parameters::Vector{T},
                                      YY::Matrix{T};
                                      mh::Bool = false,
                                      catch_errors::Bool = false)
    catch_errors = catch_errors | mh
    if mh
        try
            update!(model, parameters)
        catch err
            @printf "There was an error of type %s :" typeof(err)
            @printf "%s\n" err.msg
            return Posterior()
        end
    else
        update!(model, parameters)
    end
    return posterior(model, YY; mh=mh, catch_errors=catch_errors)

end

#=
doc"""
likelihood{T<:AbstractFloat}(model::AbstractDSGEModel, YY::Matrix{T}; mh::Bool = false, catch_errors::Bool = false)

### Parameters
-`model`: The model object
-`YY`: matrix of data for observables

## Optional Arguments
-`mh`: Whether metropolis_hastings is the caller. If `mh=true`, the transition matrices for the zero-lower-bound period are returned in a dictionary.
-`catch_errors`: If `mh = true`, `GensysErrors` should always be caught.

### Description
`likelihood` is a dsge likelihood function that can handle 2-part estimation where the observed sample contains both a normal stretch of time (in which interest rates are positive) and a stretch of time in which interest rates reach the zero lower bound. If there is a zero-lower-bound period, then we filter over the 2 periods separately. Otherwise, we filter over the main sample all at once.
"""
=#
function likelihood{T<:AbstractFloat}(model::AbstractDSGEModel{T}, YY::Matrix{T}; mh::Bool = false, catch_errors::Bool = false)
    catch_errors = catch_errors | mh
    LIKE_NULL_DICT   = Dict{Symbol, Matrix{T}}()
    LIKE_NULL_OUTPUT = (-Inf, LIKE_NULL_DICT)

    # During Metropolis-Hastings, return -∞ if any parameters are not within their bounds
    if mh
        for θ in model.parameters
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

    n_T0         = num_presample_periods(model)
    n_ant        = num_anticipated_shocks(model)
    n_ant_lags   = num_anticipated_lags(model)
    n_obs_no_ant = num_observables(model) - num_anticipated_shocks(model)
    n_obs        = num_observables(model)
    n_exo        = num_shocks_exogenous(model)

    n_states_no_ant = num_states_augmented(model) - num_anticipated_shocks(model)
    n_states_aug    = num_states_augmented(model)
    n_states        = num_states(model)
    mt_num_states   = [n_states_no_ant, n_states_no_ant, n_states_aug]

    R1[:YY] = YY[1:n_T0, 1:n_obs_no_ant]
    R2[:YY] = YY[(n_T0+1):(end-n_ant_lags-1), 1:n_obs_no_ant]
    R3[:YY] = YY[(end-n_ant_lags):end, :]


    # Step 1: solution to DSGE model - delivers transition equation for the state variables  S_t
    # transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
    # If we are in MH, then any errors coming out of gensys should be caught and a -Inf
    # posterior should be returned.
    try
        R3[:TTT], R3[:RRR], R3[:CCC] = solve(model)
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
    Measurement_R2 = measurement(model, R2[:TTT], R2[:RRR], R2[:CCC]; shocks=false)
    Measurement_R3 = measurement(model, R3[:TTT], R3[:RRR], R3[:CCC]; shocks=true)
    for d in (:ZZ, :DD, :QQ, :VVall)
        R2[d] = Measurement_R2[d]
        R3[d] = Measurement_R3[d]
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
    R1[:P0]         = dlyap(R1[:TTT], R1[:RRR]*R1[:QQ]*R1[:RRR]')
    out             = kalman_filter(R1[:YY]', 1, zeros(T, n_states_no_ant, 1), R1[:TTT], R1[:DD], R1[:ZZ], R1[:VVall], R1[:A0], R1[:P0])
    regime_likes[1] = out[:L]
    R1[:zend]       = out[:zend]
    R1[:Pend]       = out[:Pend]

    # Run Kalman filter on normal period
    zprev           = R1[:zend]
    Pprev           = R1[:Pend]
    out             = kalman_filter(R2[:YY]', 1, zeros(mt_num_states[2], 1), R2[:TTT], R2[:DD], R2[:ZZ], R2[:VVall], zprev, Pprev)
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

    out             = kalman_filter(R3[:YY]', 1, zeros(mt_num_states[3], 1), R3[:TTT], R3[:DD], R3[:ZZ], R3[:VVall], zprev, Pprev)
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



#=
doc"""
DLYAP   Discrete Lyapunov equation solver.

   X = DLYAP(A,Q) solves the discrete Lyapunov equation:

    A*X*A' - X + Q = 0

   See also  LYAP.

 This is a port of the following Matlab implementation:
 J.N. Little 2-1-86, AFP 7-28-94
 Copyright 1986-2001 The MathWorks, Inc.
 $Revision: 1.11 $  $Date: 2001/01/18 19:50:01 $

LYAP  Solve continuous-time Lyapunov equations.

  x = LYAP(a,c) solves the special form of the Lyapunov matrix
  equation:

          a*x + x*a' = -c

  See also  DLYAP.

 This is a port of the following Matlab implementation:
 S.N. Bangert 1-10-86
 Copyright 1986-2001 The MathWorks, Inc.
 $Revision: 1.10 $  $Date: 2001/01/18 19:50:23 $
 Last revised JNL 3-24-88, AFP 9-3-95

How to prove the following conversion is true.  Re: show that if
        (1) Ad X Ad' + Cd = X             Discrete lyaponuv eqn
        (2) Ac = inv(Ad + I) (Ad - I)     From dlyap
        (3) Cc = (I - Ac) Cd (I - Ac')/2  From dlyap
Then
        (4) Ac X + X Ac' + Cc = 0         Continuous lyapunov

Step 1) Substitute (2) into (3)
        Use identity 2*inv(M+I) = I - inv(M+I)*(M-I)
                                = I - (M-I)*inv(M+I) to show
        (5) Cc = 4*inv(Ad + I)*Cd*inv(Ad' + I)
Step 2) Substitute (2) and (5) into (4)
Step 3) Replace (Ad - I) with (Ad + I -2I)
        Replace (Ad' - I) with (Ad' + I -2I)
Step 4) Multiply through and simplify to get
        X -inv(Ad+I)*X -X*inv(Ad'+I) +inv(Ad+I)*Cd*inv(Ad'+I) = 0
Step 5) Left multiply by (Ad + I) and right multiply by (Ad' + I)
Step 6) Simplify to (1)
"""
=#
function dlyap(a, c)
    return dlyap!(copy(a), copy(c))
end
function dlyap!(a, c)
    m, n = size(a)
    a = (a + UniformScaling(1))\(a - UniformScaling(1))
    c = (UniformScaling(1)-a)*c*(UniformScaling(1)-a')/2

    mc, nc = size(c)

    # a and c must be square and the same size
    if (m != n) || (m != mc) || (n != nc)
        error("Dimensions do not agree.")
    elseif m == 0
        x = zeros(m, m)
        return x
    end

    # Perform schur decomposition on a (and convert to complex form)
    ta, ua, _ = schur(complex(a)) # matlab code converts to complex - come back to
    # ua, ta = rsf2csf(ua, ta)
    # Schur decomposition of a' can be calculated from that of a.
    j = m:-1:1
    ub = ua[:, j]
    tb = ta[j ,j]'

    # Check all combinations of ta(i, i)+tb(j, j) for zero
    p1 = diag(ta).' # Use .' instead of ' in case a and a' are not real
    p2 = diag(tb)
    p_sum = abs(p1) .+ abs(p2)
    if any(p_sum .== 0) || any(abs(p1 .+ p2) .< 1000*eps()*p_sum)
        error("Solution does not exist or is not unique.")
    end

    # Transform c
    ucu = -ua'*c*ub

    # Solve for first column of transformed solution
    y = complex(zeros(n, n))
    y[:, 1] = (ta + UniformScaling(tb[1, 1]))\ucu[:, 1]

    # Solve for remaining columns of transformed solution
    for k=1:n-1
        km1 = 1:k
        y[:, k+1] = (ta + UniformScaling(tb[k+1, k+1]))\(ucu[:, k+1] - y[:, km1]*tb[km1, k+1])
    end

    # Find untransformed solution
    x = ua*y*ub'

    # Ignore complex part if real inputs (better be small)
    if isreal(a) && isreal(c)
        x = real(x)
    end

    # Force x to be symmetric if c is symmetric
    if issym(c)
        x = (x+x')/2
    end

    return x
end
