#=
This code is loosely based on a routine originally copyright Federal Reserve Bank of Atlanta
and written by Iskander Karibzhanov.
=#

"""
```
kalman_filter(data, lead, CCC, TTT, DD, ZZ, VVall, z0, vz0, Ny0; allout=false)
kalman_filter(data, lead, CCC, TTT, DD, ZZ, VVall, Ny0=0; allout=false)
```

### Inputs

- `data`: a `Ny x T` matrix containing data `y(1), ... , y(T)`.
- `lead`: the number of steps to forecast after the end of the data.
- `CCC`: an `Nz x 1` vector for a time-invariant input vector in the transition equation.
- `TTT`: an `Nz x Nz` matrix for a time-invariant transition matrix in the transition
  equation.
- `DD`: an `Ny x 1` vector for a time-invariant input vector in the measurement equation.
- `ZZ`: an `Ny x Nz` matrix for a time-invariant measurement matrix in the measurement
  equation.
- `VVall`: an `Ny + Nz` x `Ny + Nz` matrix for a time-invariant variance matrix for the
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
z(t+1) = CCC+TTT*z(t)+η(t)     (state or transition equation)
y(t) = DD+ZZ*z(t)+ϵ(t)       (observation or measurement equation)
```

When `z0` and `Vz0` are omitted, the initial state vector and its covariance matrix of the
time invariant Kalman filters are computed under the stationarity condition:
```
z0 = (I-TTT)\CCC
vz0 = (I-kron(TTT,TTT))\(V(:),Nz,Nz)
```
where `TTT` and `V` are the time invariant transition matrix and the covariance matrix of
transition equation noise, and `vec(V)` is an `Nz^2` x `1` column vector that is constructed
by stacking the `Nz` columns of `V`.  Note that all eigenvalues of `TTT` are inside the unit
circle when the state space model is stationary.  When the preceding formula cannot be
applied, the initial state vector estimate is set to `a` and its covariance matrix is given
by `1E6I`.  Optionally, you can specify initial values.
"""
function kalman_filter{S<:AbstractFloat}(data::Matrix{S},
                                      TTT::Matrix{S},
                                      CCC::Vector{S},
                                      ZZ::Matrix{S},
                                      DD::Vector{S},
                                      VVall::Matrix{S},
                                      z0::Vector{S} = Vector{S}(),
                                      vz0::Matrix{S} = Matrix{S}();
                                      lead::Int = 0,
                                      allout::Bool = false,
                                      Ny0::Int = 0)
    T = size(data, 2)
    Nz = length(CCC)
    Ny = length(DD)
    V = VVall[1:Nz, 1:Nz]

    if isempty(z0) || isempty(vz0)
        e, _ = eig(TTT)
        if countnz(e*e' - eye(Nz)) == Nz^2
            z0 = (eye(Nz) - TTT)\CCC
            vz0 = reshape((eye(Nz^2)-kron(TTT,TTT))\V, Nz, Nz)
        else
            z0 = CCC
            vz0 = eye(Nz)*1e6
        end
    end

    z = z0
    P = vz0

    # Check input matrix dimensions
    @assert size(data, 1) == Ny
    @assert size(TTT) == (Nz, Nz)
    @assert size(ZZ) == (Ny, Nz)
    @assert size(VVall) == (Ny + Nz, Ny + Nz)
    @assert length(z) == Nz
    @assert size(P) == (Nz, Nz)

    # V(t) and R(t) are variances of η(t) and ϵ(t), respectively, and G(t) is a covariance
    # of η(t) and ϵ(t)
    # In dsgelh :
    # --- V is same as QQ
    # --- R is same as EE
    # --- G is same as VV = QQ*MM
    V = VVall[1:Nz, 1:Nz]
    R = VVall[(Nz+1):end, (Nz+1):end]
    G = VVall[1:Nz, (Nz+1):end]

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
        ZZ_t = ZZ[nonmissing, :]           # ZZ_t is matrix mapping states to observables
        G_t = G[:, nonmissing]             # G_t = Cov(η_t, ϵ_t)
        R_t = R[nonmissing, nonmissing]    # R_t = Var(ϵ_t)
        Ny_t = length(data_t)              # Ny_t = T is length of time
        DD_t = DD[nonmissing]                # DD_t


        ## forecasting
        z = CCC + TTT*z                    # z_{t|t-1} = CCC + TTT(Θ)*z_{t-1|t-1}
        P = TTT*P*TTT' + V                 # P_{t|t-1} = TTT(Θ)*P_{t-1|t-1}*TTT(Θ)' + TTT(Θ)*Var(η_t)*TTT(Θ)'
        dy = data_t - ZZ_t*z - DD_t        # dy = y_t - ZZ(Θ)*z_{t|t-1} - DD is prediction error or innovation
        ZG = ZZ_t*G_t                      # ZG is ZZ*Cov(η_t, ϵ_t)
        D = ZZ_t*P*ZZ_t' + ZG + ZG' + R_t  # D = ZZ*P_{t+t-1}*ZZ' + ZG + ZG' + R_t
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
        PZG = P*ZZ_t' + G_t
        z = z + PZG*ddy                    # z_{t|t} = z_{t|t-1} + P_{t|t-1}*ZZ(Θ)' + ...
        P = P - PZG/D*PZG'                 # P_{t|t} = P_{t|t-1} - PZG*(1/D)*PZG

        if allout
            PZZ = P*ZZ_t'
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
            z = TTT*z + CCC
            P = TTT*P*TTT' + V
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

"""
```
function kalman_filter_2part{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, TTT::Matrix{S} = Matrix{S}(0, 0), RRR::Matrix{S} =
    Matrix{S}(0, 0), CCC::Vector{S} = Vector{S}(0,), z0::Array{S} =
    Array{S}(0,), vz0::Matrix{S} = Matrix{S}(0, 0); DD::Array{S} = Array{S}(0,),
    lead::Int = 0, allout::Bool = false, catch_errors::Bool = false,
    include_presample::Bool = false)
```

Implements the Kalman filter, accounting for the zero lower bound.

### Inputs

- `m`: model object
- `data`: a `T x Ny` matrix containing data `y(1), ... , y(T)`.
- `TTT`: an optional `Nz x Nz` matrix for a time-invariant transition matrix in
  the transition equation. If not provided, it will be calculated.
- `RRR`: an optional `Nz` x `Nz` matrix for a time-invariant variance matrix for
  the error in the transition equation.  If not provided, it will be calculated.
- `CCC`: an `Nz` x 1` vector for a time-invariant input vector in the transition
  equation.  If not provided, it will be calculated.
- `z0`: an optional `Nz x 1` initial state vector.
- `vz0`: an optional `Nz x Nz` covariance matrix of an initial state vector.

Where:
- `Nz`: number of states
- `Ny`: number of observables
- `T`: number of periods of data

#### Keyword arguments

- `DD`: optional override for the constant term in the measurement equation. We
  use this in `durbin_koopman_smoother` because the data matrix we pass in is
  related to the states by a measurement equation with no constant term.
- `lead`: the number of steps to forecast after the end of the data.
- `allout`: indicates whether we want optional output variables returned as well
- `include_presample`: indicates whether to include presample periods in the
  returned Kalman object. If true, we concatenate Kalman objects from all three
  regimes; else only R2 and R3.

### Outputs

- a `Kalman` object. See documentation for `Kalman`.
"""
function kalman_filter_2part{S<:AbstractFloat}(m::AbstractModel,
                                               data::Matrix{S},
                                               TTT::Matrix{S} = Matrix{S}(0, 0),
                                               RRR::Matrix{S} = Matrix{S}(0, 0),
                                               CCC::Vector{S} = Vector{S}(0,),
                                               z0::Array{S}   = Array{S}(0,),
                                               vz0::Matrix{S} = Matrix{S}(0, 0);
                                               DD::Array{S}   = Array{S}(0,),
                                               lead::Int      = 0,
                                               allout::Bool   = false,
                                               catch_errors::Bool = false,
                                               include_presample::Bool = false)
    
    # Partition sample into three regimes, and store associated matrices:
    # - R1: presample
    # - R2: normal
    # - R3: zero lower bound and beyond
    R1 = Dict{Symbol, Array{S}}()
    R2 = Dict{Symbol, Array{S}}()
    R3 = Dict{Symbol, Array{S}}()

    n_T0  = n_presample_periods(m)
    n_ant = n_anticipated_shocks(m)

    n_obs_no_ant = n_observables(m) - n_anticipated_shocks(m)
    n_states_no_ant = n_states_augmented(m) - n_anticipated_shocks(m)
    n_states_aug    = n_states_augmented(m)
    nstates         = n_states(m)
    regime_states   = [n_states_no_ant, n_states_no_ant, n_states_aug]

    state_inds = inds_states_no_ant(m)
    shock_inds = inds_shocks_no_ant(m)
    obs_inds   = inds_obs_no_ant(m)

    R1[:data] = data[inds_presample_periods(m), obs_inds]
    R2[:data] = data[inds_prezlb_periods(m),    obs_inds]
    R3[:data] = data[inds_zlb_periods(m),       :]

    # Step 1: Compute the transition equation:
    #   S_t = CCC + TTT*S_{t-1} + RRR*ε_t
    # where
    #   Var(ε_t) = QQ
    # If we are in Metropolis-Hastings, then any errors coming out of `gensys`
    # should be caught and a -Inf posterior should be returned.
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
    R2[:TTT] = R3[:TTT][state_inds, state_inds]
    R2[:RRR] = R3[:RRR][state_inds, shock_inds]
    R2[:CCC] = R3[:CCC][state_inds]

    # Step 2: Define the measurement equation:
    #   Y_t = DD + ZZ*S_t + u_t
    # where
    #   Var(η_t) = EE
    #   u_t = η_t + MM*ε_t
    #   Var(u_t) = HH = EE+MM QQ MM'
    #   Cov(ε_t,u_t) = VV = QQ*MM'

    # Get measurement equation matrices set up for normal and zlb periods
    measurement_R2 = measurement(m, R2[:TTT], R2[:RRR], R2[:CCC]; shocks = false)
    measurement_R3 = measurement(m, R3[:TTT], R3[:RRR], R3[:CCC]; shocks = true)
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

    # Presample measurement & transition equation matrices are same as normal
    # period
    for d in (:TTT, :RRR, :QQ, :ZZ, :DD, :VVall)
        R1[d] = R2[d]
    end

    # Step 3: Compute log-likelihood using the Kalman filter.
    # Note that `kalman_filter` assumes a transition equation of the form:
    #   S_t = TTT*S_{t-1} + ε2_t
    # where ε2_t = RRR*ε_t. Therefore redefine:
    #   QQ2 = Var(ε2_t) = RRR*QQ*RRR'
    #   VV2 = Cov(ε2_t, u_t) = RRR*VV
    #   VVall = Var([ε2_t; u_t])    (joint variance of the two shocks)

    # Run Kalman filter on presample
    R1[:A0] = if isempty(z0)
        zeros(S, n_states_no_ant)
    else
        z0[state_inds]
    end
    R1[:P0] = solve_discrete_lyapunov(R1[:TTT], R1[:RRR]*R1[:QQ]*R1[:RRR]')
    k1 = kalman_filter(R1[:data]', R1[:TTT], zeros(S, regime_states[1]),
        R1[:ZZ], R1[:DD], R1[:VVall], R1[:A0], R1[:P0]; lead = 1, allout = allout)

    # Run Kalman filter on normal period
    k2 = kalman_filter(R2[:data]', R2[:TTT], zeros(regime_states[2]),
        R2[:ZZ], R2[:DD], R2[:VVall], k1[:zend], k1[:Pend]; lead = 1, allout = allout)

    # Run Kalman filter on ZLB period
    zprev = zeros(S, n_states_aug)
    Pprev = zeros(S, n_states_aug, n_states_aug)
    zprev[state_inds] = k2[:zend]
    Pprev[state_inds, state_inds] = k2[:Pend]
    k3 = kalman_filter(R3[:data]', R3[:TTT], zeros(regime_states[3]),
        R3[:ZZ], R3[:DD], R3[:VVall], zprev, Pprev; lead = 1, allout = allout)

    # Concatenate Kalman objects
    if include_presample
        k12 = cat(m, k1, k2; allout = allout)
        k = cat(m, k12, k3, regime_switch = true, allout = allout)
    else
        k = cat(m, k2, k3; regime_switch = true, allout = allout)
    end

    ## Return concatenated Kalman and system matrices for each regime
    return k, R1, R2, R3
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
    vpred::Array{S, 3}    # predicted variance-covariance matrix for states in period T+1
    yprederror::Matrix{S} 
    ystdprederror::Matrix{S}
    rmse::Matrix{S}
    rmsd::Matrix{S}        # 
    filt::Matrix{S}        # filtered states
    vfilt::Array{S, 3}     # mean square errors of filtered state vectors
    z0::Vector{S}          # starting-period state vector
    vz0::Matrix{S}         # starting-period variance-covariance matrix for the states
end

function Kalman{S<:AbstractFloat}(L::S,
                                  zend::Vector{S},
                                  Pend::Matrix{S},
                                  pred::Matrix{S}          = Matrix{S}(0, 0),
                                  vpred::Array{S, 3}       = Array{S}(0, 0, 0),
                                  yprederror::Matrix{S}    = Matrix{S}(0, 0),
                                  ystdprederror::Matrix{S} = Matrix{S}(0, 0),
                                  rmse::Matrix{S}          = Matrix{S}(0, 0),
                                  rmsd::Matrix{S}          = Matrix{S}(0, 0),
                                  filt::Matrix{S}          = Matrix{S}(0, 0),
                                  vfilt::Array{S, 3}       = Array{S}(0, 0, 0),
                                  z0::Vector{S}            = Vector{S}(0),
                                  vz0::Matrix{S}           = Matrix{S}(0, 0))
    return Kalman{S}(L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt, z0, vz0)
end

function Base.getindex(K::Kalman, d::Symbol)
    if d in (:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse, :rmsd,
             :filt, :vfilt, :z0, :vz0)
        return getfield(K, d)
    else
        throw(KeyError(d))
    end
end

function Base.cat{S<:AbstractFloat}(m::AbstractModel, k1::Kalman{S},
    k2::Kalman{S}; regime_switch::Bool = false, allout::Bool = false)

    # If k1 results from calling the Kalman filter on pre-ZLB data and k2 from
    # calling it on data under the ZLB, then we must augment the fields of k1 to
    # accommodate the states and observables corresponding to anticipated policy
    # shocks

    # However, if !allout, then we don't need to augment anything because we'll
    # only fill out L, zend, and Pend
    if regime_switch && allout

        n_states_aug = n_states_augmented(m)
        n_obs        = n_observables(m)
        n_k1_periods = size(k1[:pred], 2)

        # Initialize fields for augmented k1
        pred = zeros(S, n_states_aug, n_k1_periods)
        vpred = zeros(S, n_states_aug, n_states_aug, n_k1_periods)
        yprederror = zeros(S, n_obs, n_k1_periods)
        ystdprederror = zeros(S, n_obs, n_k1_periods)
        rmse = zeros(S, 1, n_obs)
        rmsd = zeros(S, 1, n_obs)
        filt = zeros(S, n_states_aug, n_k1_periods)
        vfilt = zeros(S, n_states_aug, n_states_aug, n_k1_periods)
        z0 = zeros(S, n_states_aug)
        vz0 = zeros(S, n_states_aug, n_states_aug)
        k1_new = Kalman(k1[:L], k1[:zend], k1[:Pend], pred, vpred, yprederror,
            ystdprederror, rmse, rmsd, filt, vfilt, z0, vz0)

        state_inds = inds_states_no_ant(m)
        obs_inds   = inds_obs_no_ant(m)

        # Augment fields
        k1_new[:pred][state_inds, :] = k1[:pred]
        k1_new[:vpred][state_inds, state_inds, :] = k1[:vpred]
        k1_new[:yprederror][obs_inds, :] = k1[:yprederror]
        k1_new[:ystdprederror][obs_inds, :] = k1[:ystdprederror]
        k1_new[:rmse][:, obs_inds] = k1[:rmse]
        k1_new[:rmsd][:, obs_inds] = k1[:rmsd]
        k1_new[:filt][state_inds, :] = k1[:pred]
        k1_new[:vfilt][state_inds, state_inds, :] = k1[:vfilt]
        k1_new[:z0][state_inds] = k1[:z0]
        k1_new[:vz0][state_inds, state_inds] = k1[:vz0]

        # Replace k1 with augmented fields
        k1 = k1_new
    end

    L = k1[:L] + k2[:L]
    zend = k2[:zend]
    Pend = k2[:Pend]

    if allout
        pred = hcat(k1[:pred], k2[:pred])
        vpred = cat(3, k1[:vpred], k2[:vpred])
        yprederror = hcat(k1[:yprederror], k2[:yprederror])
        ystdprederror = hcat(k1[:ystdprederror], k2[:yprederror])
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))
        filt = hcat(k1[:filt], k2[:filt])
        vfilt = cat(3, k1[:vfilt], k2[:vfilt])
        z0   = k1[:z0]
        vz0  = k1[:vz0]

        return Kalman(L, zend, Pend, pred, vpred, yprederror, ystdprederror,
            rmse, rmsd, filt, vfilt, z0, vz0)
    else
        return Kalman(L, zend, Pend)
    end
end