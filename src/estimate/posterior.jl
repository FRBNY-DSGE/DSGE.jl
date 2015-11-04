using Compat

#=
doc"""
Calculates (log of) the joint density of the model parameters.
"""
=#
function prior(model::AbstractDSGEModel)
    x = zero(Float64)
    for θ in model.parameters
        if !θ.fixed
            x += logpdf(θ)
        end
    end
    return x
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
function posterior{T<:AbstractFloat}(model::AbstractDSGEModel, YY::Matrix{T}; mh::Bool = false, catch_errors::Bool = false)
    if mh
        catch_errors = true
        like, out = likelihood(model, YY; mh=mh)
        post = like + prior(model)
        return post, like, out
    else
        return likelihood(model, YY; mh=mh, catch_errors=catch_errors) + prior(model)
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
function posterior!{T<:AbstractFloat}(model::AbstractDSGEModel, parameters::Vector{T}, YY::Matrix{T}; mh::Bool = false, catch_errors::Bool = false)
    MH_NULL_OUTPUT = (-Inf, Dict{Symbol, Any}())
        
    if mh

        try
            update!(model, parameters)
        catch err
            @printf "There was an error of type %s :" typeof(err)
            @printf "%s\n" err.msg
            return MH_NULL_OUTPUT
        end

        return posterior(model, YY; mh=true, catch_errors=true)
    else
        return posterior(model, YY; mh=mh, catch_errors=catch_errors)
    end
    

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
function likelihood{T<:AbstractFloat}(model::AbstractDSGEModel, YY::Matrix{T}; mh::Bool = false, catch_errors::Bool = false)
    MH_NULL_OUTPUT = (-Inf, Dict{Symbol, Any}())
    GENSYS_ERROR_OUTPUT = -Inf

    # During Metropolis-Hastings, return -∞ if any parameters are not within their bounds
    if mh
        catch_errors = true  # for consistency
        for θ in model.parameters
            (left, right) = θ.valuebounds
            if !θ.fixed && !(left <= θ.value <= right)
                return MH_NULL_OUTPUT
            end
        end
    end

    # Partition sample into presample, normal, and zero lower bound periods
    presample = Dict{Symbol, Any}()
    normal = Dict{Symbol, Any}()
    zlb = Dict{Symbol, Any}()
    mt = [presample, normal, zlb]

    presample[:num_observables] = num_observables(model) - num_anticipated_shocks(model)
    normal[:num_observables] = num_observables(model) - num_anticipated_shocks(model)
    zlb[:num_observables] = num_observables(model)

    presample[:num_states] = num_states_augmented(model) - num_anticipated_shocks(model)
    normal[:num_states] = num_states_augmented(model) - num_anticipated_shocks(model)
    zlb[:num_states] = num_states_augmented(model)

    presample[:YY] = YY[1:num_presample_periods(model), 1:presample[:num_observables]]
    normal[:YY] = YY[(num_presample_periods(model)+1):(end-num_anticipated_lags(model)-1), 1:normal[:num_observables]]
    zlb[:YY] = YY[(end-num_anticipated_lags(model)):end, :]


    # Step 1: solution to DSGE model - delivers transition equation for the state variables  S_t
    # transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
    # If we are in MH, then any errors coming out of gensys should be caught and a -Inf
    # posterior should be returned.
    try
        zlb[:TTT], zlb[:RRR], zlb[:CCC] = solve(model)
    catch err
        if catch_errors && isa(err, GensysError)
            info(err.msg)
            if mh
                return MH_NULL_OUTPUT
            else
                return GENSYS_ERROR_OUTPUT
            end
        else
            rethrow(err)
        end
    end

    # Get normal, no ZLB m matrices
    state_inds = [1:(num_states(model)-num_anticipated_shocks(model)); (num_states(model)+1):num_states_augmented(model)]
    shock_inds = 1:(num_shocks_exogenous(model)-num_anticipated_shocks(model))

    normal[:TTT] = zlb[:TTT][state_inds, state_inds]
    normal[:RRR] = zlb[:RRR][state_inds, shock_inds]
    normal[:CCC] = zlb[:CCC][state_inds, :]



    ## step 2: define the measurement equation: X_t = ZZ*S_t + D + u_t
    ## where u_t = eta_t+MM* eps_t with var(eta_t) = EE
    ## where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

    # Get measurement equation matrices set up for all periods that aren't the presample
    for p = 2:3
        #try
            shocks = (p == 3)
            mt[p][:ZZ], mt[p][:DD], mt[p][:QQ], mt[p][:EE], mt[p][:MM] = measurement(model, mt[p][:TTT], mt[p][:RRR], mt[p][:CCC]; shocks=shocks)
        #catch
            # Error thrown during gensys
            #    return -Inf
        #end

        mt[p][:HH] = mt[p][:EE] + mt[p][:MM]*mt[p][:QQ]*mt[p][:MM]'
        mt[p][:VV] = mt[p][:QQ]*mt[p][:MM]'
        mt[p][:VVall] = [[mt[p][:RRR]*mt[p][:QQ]*mt[p][:RRR]'   mt[p][:RRR]*mt[p][:VV]];
                          [mt[p][:VV]'*mt[p][:RRR]'               mt[p][:HH]]]
    end

    # TODO: Incorporate this into measurement equation (why is this only done for normal period?)
    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(normal[:CCC] != 0)
        normal[:DD] += normal[:ZZ]*((UniformScaling(1) - normal[:TTT])\normal[:CCC])
    end

    # Presample measurement & transition equation matrices are same as normal period
    presample[:TTT], presample[:RRR], presample[:QQ] = normal[:TTT], normal[:RRR], normal[:QQ]
    presample[:ZZ], presample[:DD], presample[:VVall] = normal[:ZZ], normal[:DD], normal[:VVall]



    ## step 3: compute log-likelihood using Kalman filter - written by Iskander
    ##         note that Iskander's program assumes a transition equation written as:
    ##         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t
    ##         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
    ##         and  VV2 = cov(eps2_t,u_u) = RRR*VV
    ##         define VVall as the joint variance of the two shocks VVall = var([eps2_tu_t])

    # Run Kalman filter on presample
    presample[:A0] = zeros(presample[:num_states], 1)
    presample[:P0] = dlyap!(copy(presample[:TTT]), copy(presample[:RRR]*presample[:QQ]*presample[:RRR]'))
    presample[:pyt], presample[:zend], presample[:Pend] = kalcvf2NaN(presample[:YY]', 1, zeros(presample[:num_states], 1), presample[:TTT], presample[:DD], presample[:ZZ], presample[:VVall], presample[:A0], presample[:P0])

    # Run Kalman filter on normal and ZLB periods
    for p = 2:3
        if p == 2
            zprev = presample[:zend]
            Pprev = presample[:Pend]
        else
            # This section expands the number of states to accomodate extra states for the
            # anticipated policy shocks. It does so by taking the zend and Pend for the
            # state space without anticipated policy shocks, then shoves in nant
            # zeros in the middle of zend and Pend in the location of
            # the anticipated shock entries.
            before_shocks = 1:(num_states(model)-num_anticipated_shocks(model))
            after_shocks_old = (num_states(model)-num_anticipated_shocks(model)+1):(num_states_augmented(model)-num_anticipated_shocks(model))
            after_shocks_new = (num_states(model)+1):num_states_augmented(model)

            zprev = [normal[:zend][before_shocks, :];
                     zeros(num_anticipated_shocks(model), 1);
                     normal[:zend][after_shocks_old, :]]

            Pprev = zeros(num_states_augmented(model), num_states_augmented(model))
            Pprev[before_shocks, before_shocks] = normal[:Pend][before_shocks, before_shocks]
            Pprev[before_shocks, after_shocks_new] = normal[:Pend][before_shocks, after_shocks_old]
            Pprev[after_shocks_new, before_shocks] = normal[:Pend][after_shocks_old, before_shocks]
            Pprev[after_shocks_new, after_shocks_new] = normal[:Pend][after_shocks_old, after_shocks_old]
        end

        mt[p][:pyt], mt[p][:zend], mt[p][:Pend] = kalcvf2NaN(mt[p][:YY]', 1, zeros(mt[p][:num_states], 1), mt[p][:TTT], mt[p][:DD], mt[p][:ZZ], mt[p][:VVall], zprev, Pprev)
    end

    # Return total log-likelihood, excluding the presample
    like = normal[:pyt] + zlb[:pyt]
    if mh
        return like, zlb
    else
        return like
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
