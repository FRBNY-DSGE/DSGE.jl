#=
This code is loosely based on a routine originally copyright Federal Reserve Bank of Atlanta
and written by Iskander Karibzhanov.
=#

"""
```
kalman_filter(data, lead, a, F, b, H, var, z0, vz0, Ny0; allout=false)
kalman_filter(data, lead, a, F, b, H, var, Ny0=0; allout=false)
```

### Inputs

- `data`: a `Ny x T` matrix containing data `y(1), ... , y(T)`.
- `lead`: the number of steps to forecast after the end of the data.
- `a`: an `Nz x 1` vector for a time-invariant input vector in the transition equation.
- `F`: an `Nz x Nz` matrix for a time-invariant transition matrix in the transition
  equation.
- `b`: an `Ny x 1` vector for a time-invariant input vector in the measurement equation.
- `H`: an `Ny x Nz` matrix for a time-invariant measurement matrix in the measurement
  equation.
- `var`: an `Ny + Nz` x `Ny + Nz` matrix for a time-invariant variance matrix for the
  error in the transition equation and the error in the measurement equation, that is,
  `[η(t)', ϵ(t)']'`.

#### Optional Inputs
- `z0`: an optional `Nz x 1` initial state vector.
- `vz0`: an optional `Nz x Nz` covariance matrix of an initial state vector.
- `Ny0`: an optional scalar indicating the number of periods of presample
  (i.e. the number of periods which we don't add to the likelihood). If `Ny0 >
  0`, then we also set the `z0` and `vz0` fields in the returned `Kalman` object
  to be the states and variance-covariance matrices at the end of the
  presample/beginning of the main sample
- `allout`: an optional keyword argument indicating whether we want optional output
  variables returned as well

Where:
- `Nz`: number of states
- `Ny`: number of observables
- `T`: number of time periods for which we have data

### Outputs

- a `Kalman` object. See documentation for `Kalman`.


### Notes

The state space model is defined as follows:
```
z(t+1) = a+F*z(t)+η(t)     (state or transition equation)
y(t) = b+H*z(t)+ϵ(t)       (observation or measurement equation)
```

When `z0` and `Vz0` are omitted, the initial state vector and its covariance matrix of the
time invariant Kalman filters are computed under the stationarity condition:
```
z0 = (I-F)\a
vz0 = (I-kron(F,F))\(V(:),Nz,Nz)
```
where `F` and `V` are the time invariant transition matrix and the covariance matrix of
transition equation noise, and `vec(V)` is an `Nz^2` x `1` column vector that is constructed
by stacking the `Nz` columns of `V`.  Note that all eigenvalues of `F` are inside the unit
circle when the state space model is stationary.  When the preceding formula cannot be
applied, the initial state vector estimate is set to `a` and its covariance matrix is given
by `1E6I`.  Optionally, you can specify initial values.
"""
function kalman_filter{S<:AbstractFloat}(data::Matrix{S},
                                      lead::Int64,
                                      a::Vector{S},
                                      F::Matrix{S},
                                      b::Vector{S},
                                      H::Matrix{S},
                                      var::Matrix{S},
                                      z0::Vector{S},
                                      vz0::Matrix{S},
                                      Ny0::Int = 0;
                                      allout::Bool = false)
    T = size(data, 2)
    Nz = length(a)
    Ny = length(b)

    z = z0
    P = vz0

    # Check input matrix dimensions
    @assert size(data, 1) == Ny
    @assert size(F) == (Nz, Nz)
    @assert size(H) == (Ny, Nz)
    @assert size(var) == (Ny + Nz, Ny + Nz)
    @assert length(z) == Nz
    @assert size(P) == (Nz, Nz)

    # V(t) and R(t) are variances of η(t) and ϵ(t), respectively, and G(t) is a covariance
    # of η(t) and ϵ(t)
    # In dsgelh :
    # --- V is same as QQ
    # --- R is same as EE
    # --- G is same as VV = QQ*MM
    V = var[1:Nz, 1:Nz]
    R = var[(Nz+1):end, (Nz+1):end]
    G = var[1:Nz, (Nz+1):end]

    if allout
        pred          = zeros(S, Nz, T)
        vpred         = zeros(S, Nz, Nz, T)
        yprederror    = NaN*zeros(S, Ny, T)
        ystdprederror = NaN*zeros(S, Ny, T)
        filt          = zeros(Nz, T)
        vfilt         = zeros(Nz, Nz, T)
    end

    L = zero(S)

    for t = 1:T
        # If an element of the vector y(t) is missing (NaN) for the observation t, the
        #   corresponding row is ditched from the measurement equation.
        nonmissing = !isnan(data[:, t])
        data_t = data[nonmissing, t]       # data_t = Y_T = [y1, y2, ..., yT] is matrix of observable data time-series
        H_t = H[nonmissing, :]             # H_t = DD is matrix mapping states to observables
        G_t = G[:, nonmissing]             # G_t = Cov(η_t, ϵ_t)
        R_t = R[nonmissing, nonmissing]    # R_t = Var(ϵ_t)
        Ny_t = length(data_t)              # Ny_t = T is length of time
        b_t = b[nonmissing]                # b_t = DD


        ## forecasting
        z = a + F*z                        # z_{t|t-1} = a + F(Θ)*z_{t-1|t-1}
        P = F*P*F' + V                     # P_{t|t-1} = F(Θ)*P_{t-1|t-1}*F(Θ)' + F(Θ)*Var(η_t)*F(Θ)'
        dy = data_t - H_t*z - b_t          # dy = y_t - H(Θ)*z_{t|t-1} - DD is prediction error or innovation
        HG = H_t*G_t                       # HG is ZZ*Cov(η_t, ϵ_t)
        D = H_t*P*H_t' + HG + HG' + R_t    # D = ZZ*P_{t+t-1}*ZZ' + HG + HG' + R_t
        D = (D+D')/2

        if allout
            pred[:, t]                   = z
            vpred[:, :, t]               = P
            yprederror[nonmissing, t]    = dy
            ystdprederror[nonmissing, t] = dy./sqrt(diag(D))
        end

        ddy = D\dy
        # We evaluate the log likelihood function by adding values of L at every iteration
        #   step (for each t = 1,2,...T)
        if t > Ny0
            L += -log(det(D))/2 - first(dy'*ddy/2) - Ny_t*log(2*pi)/2
        end

        ## updating
        PHG = P*H_t' + G_t
        z = z + PHG*ddy                    # z_{t|t} = z_{t|t-1} + P_{t|t-1}*H(Θ)' + ...
        P = P - PHG/D*PHG'                 # P_{t|t} = P_{t|t-1} - PHG*(1/D)*PHG

        if allout
            PH = P*H_t'
            filt[:, t]     = z
            vfilt[:, :, t] = P
        end

        # If Ny0 > 0 (positive number of presample periods), then we reassign
        # `z0` and `P0` to be their values at the end of the presample/beginning
        # of the main sample
        if t == Ny0
            z0 = z
            P0 = P
        end
    end

    zend = z
    Pend = P

    if allout && lead > 1
        for t = (T+2):(T+lead)
            z = F*z + a
            P = F*P*F' + V
            pred[:, t]     = z
            vpred[:, :, t] = P
        end
    end

    if allout
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))

        return Kalman(L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt, z0, vz0)
    else
        return Kalman(L, zend, Pend)
    end

end

function kalman_filter{S<:AbstractFloat}(data::Matrix{S},
                                      lead::Int64,
                                      a::Vector{S},
                                      F::Matrix{S},
                                      b::Vector{S},
                                      H::Matrix{S},
                                      var::Matrix{S},
                                      Ny0::Int = 0;
                                      allout::Bool = false)
    Nz = length(a)
    V = var[1:Nz, 1:Nz]

    e, _ = eig(F)
    if countnz(e*e' - eye(Nz)) == Nz^2
        z0 = (eye(Nz) - F)\a
        vz0 = reshape((eye(Nz^2)-kron(F,F))\V, Nz, Nz)
    else
        z0 = a
        vz0 = eye(Nz)*1e6
    end

    return kalman_filter(data, lead, a, F, b, H, var, z0, vz0, Ny0; allout=allout)
end


"""
```
kalman_filter_2part{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
     TTT::Matrix{S}, RRR::Matrix{S}, CCC::Matrix{S}, z0::Array{S},
     vz0::Matrix{S}; lead::Int, Ny0::Int, allout::Bool,
     catch_errors::Bool)
```

Implments the Kalman Filter over 2 regimes. The first regime does not
incorporate interest rate expectations data; the second does.

### Inputs
Where:
- `Nz`: number of states
- `Ny`: number of observables

- `m`: model object
- `data`: a `T x Ny` matrix containing data `y(1), ... , y(T)`.
- `TTT`: an optional `Nz x Nz` matrix for a time-invariant transition matrix in the transition
  equation. If not provided, it will be calculated.
- `RRR`: an optional `Nz` x `Nz` matrix for a time-invariant variance matrix
  for the error in the transition equation.  If not provided, it will be calculated.
- `CCC`: an `Nz x 1` vector for a time-invariant input vector in the
  transition equation.  If not provided, it will be calculated.
- `z0`: an optional `Nz x 1` initial state vector.
- `vz0`: an optional `Nz x Nz` covariance matrix of an initial state vector.

#### Keyword arguments

- `lead`: the number of steps to forecast after the end of the data.
- `Ny0`: an optional scalar indicating the number of periods of presample (i.e. the number
  of periods which we don't add to the likelihood)
- `allout`: an optional keyword argument indicating whether we want optional output
  variables returned as well
- `augment_states`: an optional keyword argument indicating whether we want to
  augment the pre-ZLB `filt`, `pred`, and `vpred` matrices with rows of zeros
  for the states corresponding to anticipated policy shocks
"""
function kalman_filter_2part{S<:AbstractFloat}(m::AbstractModel,
                                               data::Matrix{S},
                                               TTT::Matrix{S} = Matrix{S}(0, 0),
                                               RRR::Matrix{S} = Matrix{S}(0, 0),
                                               CCC::Matrix{S} = Matrix{S}(0, 0),
                                               z0::Array{S}   = Array{S}(0),
                                               vz0::Matrix{S} = Matrix{S}(0, 0);
                                               DD::Array{S}   = Array{S}(0),
                                               lead::Int      = 0,
                                               Ny0::Int       = 0,
                                               allout::Bool   = false,
                                               catch_errors::Bool = false,
                                               augment_states::Bool = false)
    
    # Partition sample into three regimes, and store associated matrices:
    # - R1: presample
    # - R2: normal
    # - R3: zero lower bound and beyond
    R1 = Dict{Symbol, Array{S}}()
    R2 = Dict{Symbol, Array{S}}()
    R3 = Dict{Symbol, Array{S}}()
    regime_mats = [R1, R2, R3]

    n_T0         = n_presample_periods(m)
    n_ant        = n_anticipated_shocks(m)

    n_obs_no_ant = n_observables(m) - n_anticipated_shocks(m)
    n_obs        = n_observables(m)
    n_exo        = n_shocks_exogenous(m)

    n_states_no_ant = n_states_augmented(m) - n_anticipated_shocks(m)
    n_states_aug    = n_states_augmented(m)
    nstates         = n_states(m)
    regime_states   = [n_states_no_ant, n_states_no_ant, n_states_aug]

    R1[:data] = data[1:(index_prezlb_start(m)-1),                   1:n_obs_no_ant]
    R2[:data] = data[index_prezlb_start(m):(index_zlb_start(m)-1), 1:n_obs_no_ant]
    R3[:data] = data[index_zlb_start(m):end,                        :]

    # Step 1: solution to DSGE model - delivers transition equation for the state variables
    # transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
    # If we are in Metropolis-Hastings, then any errors coming out of gensys should be caught and a -Inf
    # posterior should be returned.
    if isempty(TTT) || isempty(RRR) || isempty(CCC)
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
    else
        R3[:TTT], R3[:RRR], R3[:CCC] = TTT, RRR, CCC
    end

    
    # Get normal, no ZLB matrices
    state_inds = [1:(nstates-n_ant); (nstates+1):n_states_aug]
    shock_inds = 1:(n_exo-n_ant)
    obs_inds   = 1:(n_obs-n_ant)

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

    # If we pass in DD as a kwarg (most often a vector of zeros, as in the
    # Durbin-Koopman smoother), we want to use that DD instead of the one
    # calculated from the measurement equation
    if !isempty(DD)
        R2[:DD] = DD[obs_inds]
        R3[:DD] = DD
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
    R1[:A0] = if isempty(z0)
        zeros(S, n_states_no_ant)
    else
        z0[state_inds]
    end

    R1[:P0]         = solve_discrete_lyapunov(R1[:TTT], R1[:RRR]*R1[:QQ]*R1[:RRR]')

    out             = kalman_filter(R1[:data]', 1, zeros(S, regime_states[1]), R1[:TTT], R1[:DD], R1[:ZZ], R1[:VVall], R1[:A0], R1[:P0], allout=allout)

    R1[:like]       = Matrix{S}(1,1)
    R1[:like][1,1]  = out[:L]
    for d in [:zend, :Pend, :filt, :pred, :vpred]
        R1[d] = out[d]
    end

    # Run Kalman filter on normal period
    zprev           = R1[:zend]
    Pprev           = R1[:Pend]
    out             = kalman_filter(R2[:data]', 1, zeros(regime_states[2]), R2[:TTT], R2[:DD], R2[:ZZ], R2[:VVall], zprev, Pprev, allout=allout)
        
    R2[:like]       = Matrix{S}(1,1)
    R2[:like][1,1]  = out[:L]
    for d in [:zend, :Pend, :filt, :pred, :vpred, :z0, :vz0]
        R2[d] = out[d]
    end

    # Run Kalman filter on ZLB period
    # This section expands the number of states to accomodate extra states for the
    # anticipated policy shocks. It does so by taking the zend and Pend for the
    # state space without anticipated policy shocks, then shoves in nant
    # zeros in the middle of zend and Pend in the location of
    # the anticipated shock entries.
    before_shocks    = 1:(nstates-n_ant)
    after_shocks_old = (nstates-n_ant+1):(n_states_aug-n_ant)
    after_shocks_new = (nstates+1):n_states_aug

    zprev = [R2[:zend][before_shocks];
             zeros(S, n_ant);
             R2[:zend][after_shocks_old]]

    Pprev                                     = zeros(S, n_states_aug, n_states_aug)
    Pprev[before_shocks, before_shocks]       = R2[:Pend][before_shocks, before_shocks]
    Pprev[before_shocks, after_shocks_new]    = R2[:Pend][before_shocks, after_shocks_old]
    Pprev[after_shocks_new, before_shocks]    = R2[:Pend][after_shocks_old, before_shocks]
    Pprev[after_shocks_new, after_shocks_new] = R2[:Pend][after_shocks_old, after_shocks_old]

    out             = kalman_filter(R3[:data]', 1, zeros(regime_states[3]), R3[:TTT], R3[:DD], R3[:ZZ], R3[:VVall], zprev, Pprev, allout=allout)

    R3[:like]       = Matrix{S}(1,1)
    R3[:like][1,1]  = out[:L]
    for d in [:zend, :Pend, :filt, :pred, :vpred, :z0, :vz0]
        R3[d] = out[d]
    end

    # TODO: Is there some way to iterate through `R1` and `R2`?
    # If `augment_states`, then we expand the number of states in the `filt`,
    # `pred`, and `vpred` matrices before returning
    if augment_states
        R1_filt_small  = R1[:filt]
        R1[:filt]      = zeros(n_states_aug, n_T0)
        R1[:filt][before_shocks,    :] = R1_filt_small[before_shocks,    :]
        R1[:filt][after_shocks_new, :] = R1_filt_small[after_shocks_old, :]
        
        R1_pred_small  = R1[:pred]
        R1[:pred]      = zeros(n_states_aug, n_T0)
        R1[:pred][before_shocks,    :] = R1_pred_small[before_shocks,    :]
        R1[:pred][after_shocks_new, :] = R1_pred_small[after_shocks_old, :]

        R1_vpred_small = R1[:vpred]
        R1[:vpred]     = zeros(n_states_aug, n_states_aug, n_T0)
        R1[:vpred][before_shocks,    before_shocks,    :] = R1_vpred_small[before_shocks,    before_shocks,    :]
        R1[:vpred][before_shocks,    after_shocks_new, :] = R1_vpred_small[before_shocks,    after_shocks_old, :]
        R1[:vpred][after_shocks_new, before_shocks,    :] = R1_vpred_small[after_shocks_old, before_shocks,    :]
        R1[:vpred][after_shocks_new, after_shocks_new, :] = R1_vpred_small[after_shocks_old, after_shocks_old, :]

        R2_filt_small  = R2[:filt]
        R2[:filt]      = zeros(n_states_aug, n_prezlb_periods(m))
        R2[:filt][before_shocks,    :] = R2_filt_small[before_shocks,    :]
        R2[:filt][after_shocks_new, :] = R2_filt_small[after_shocks_old, :]
        
        R2_pred_small  = R2[:pred]
        R2[:pred]      = zeros(n_states_aug, n_prezlb_periods(m))
        R2[:pred][before_shocks,    :] = R2_pred_small[before_shocks,    :]
        R2[:pred][after_shocks_new, :] = R2_pred_small[after_shocks_old, :]

        R2_vpred_small = R2[:vpred]
        R2[:vpred]     = zeros(n_states_aug, n_states_aug, n_prezlb_periods(m))
        R2[:vpred][before_shocks,    before_shocks,    :] = R2_vpred_small[before_shocks,    before_shocks,    :]
        R2[:vpred][before_shocks,    after_shocks_new, :] = R2_vpred_small[before_shocks,    after_shocks_old, :]
        R2[:vpred][after_shocks_new, before_shocks,    :] = R2_vpred_small[after_shocks_old, before_shocks,    :]
        R2[:vpred][after_shocks_new, after_shocks_new, :] = R2_vpred_small[after_shocks_old, after_shocks_old, :]
end

    ## Return outputs from both regimes
    return R2, R3, R1
   
end


"""
```
Kalman{S<:AbstractFloat}
```

### Fields:
- `L`: value of the average log likelihood function of the SSM under assumption that
  observation noise ϵ(t) is normally distributed
- `zend`: state vector in the last period for which data is provided
- `Pend`: variance-covariance matrix for `zend`
- `z0`: starting-period state vector. If there are presample periods in the
  data, then `z0` is the state vector at the end of the presample/beginning of
  the main sample
- `vz0`: variance-covariance matrix for `z0`

#### Fields filled in when `allout=true` in a call to `kalman_filter`:
- `pred`: a `Nz` x `T+lead` matrix containing one-step predicted state vectors
- `vpred`: a `Nz` x `Nz` x `T+lead` matrix containing mean square errors of predicted
  state vectors

- `filt`: an `Nz` x `T` matrix containing filtered state vectors
- `vfilt`: an `Nz` x `Nz` x `T` matrix containing mean square errors of filtered
  state vectors
"""
immutable Kalman{S<:AbstractFloat}
    L::S                  # Likelihood
    zend::Vector{S}       # last-period state vector
    Pend::Matrix{S}       # last-period variance-covariance matrix for the states
    pred::Matrix{S}       # predicted value of states in period T+1
    vpred::Array{S,3}     # predicted variance-covariance matrix for states in period T+1
    yprederror::Matrix{S} 
    ystdprederror::Matrix{S}
    rmse::Matrix{S}
    rmsd::Matrix{S}        # 
    filt::Matrix{S}        # filtered states
    vfilt::Array{S,3}      # mean square errors of filtered state vectors
    z0::Vector{S}           # starting-period state vector
    vz0::Matrix{S}         # starting-period variance-covariance matrix for the states
end
function Kalman{S<:AbstractFloat}(L::S,
                                  zend::Vector{S},
                                  Pend::Matrix{S},
                                  pred::Matrix{S}          = Matrix{S}(),
                                  vpred::Array{S,3}        = Array{S}(0,0,0),
                                  yprederror::Matrix{S}    = Matrix{S}(),
                                  ystdprederror::Matrix{S} = Matrix{S}(),
                                  rmse::Matrix{S}          = Matrix{S}(),
                                  rmsd::Matrix{S}          = Matrix{S}(),
                                  filt::Matrix{S}          = Matrix{S}(),
                                  vfilt::Array{S,3}        = Array{S}(0,0,0),
                                  z0::Vector{S}            = Vector{S}(0),
                                  vz0::Matrix{S}           = Array{S}(0,0))
    return Kalman{S}(L,zend,Pend,pred,vpred,yprederror,ystdprederror,rmse,rmsd,filt,vfilt,z0,vz0)
end
function Base.getindex(K::Kalman, d::Symbol)
    if d in (:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse, :rmsd,
             :filt, :vfilt, :z0, :vz0)
        return getfield(K, d)
    else
        throw(KeyError(d))
    end
end
