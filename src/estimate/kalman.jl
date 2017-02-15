#=
This code is loosely based on a routine originally copyright Federal Reserve Bank of Atlanta
and written by Iskander Karibzhanov.
=#

"""
```
kalman_filter{S<:AbstractFloat}(data::Matrix{S},
    TTT::Matrix{S}, CCC::Vector{S}, ZZ::Matrix{S}, DD::Vector{S}, VVall::Matrix{S},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = false)
```

### Inputs

- `data`: a `Ny` x `T` `Matrix` containing data `y(1), ... , y(T)`.
- `TTT`: an `Nz` x `Nz` `Matrix` for a time-invariant transition matrix in the transition
  equation.
- `CCC`: an `Nz` x 1 `Vector` for a time-invariant input vector in the transition equation.
- `ZZ`: an `Ny` x `Nz` `Matrix` for a time-invariant measurement matrix in the measurement
  equation.
- `DD`: an `Ny` x 1 constant vector in measurement equation
- `VVall`: See `Measurement` type for description

#### Optional Inputs
- `z0`: an optional `Nz` x 1 initial state vector.
- `P0`: an optional `Nz` x `Nz` covariance matrix of an initial state vector.
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
z(t+1) = CCC + TTT*z(t) + RRR*ϵ(t)   (state or transition equation)
y(t) = DD + ZZ*z(t) + u(t)           (observation or measurement equation)
```

When `z0` and `P0` are omitted, the initial state vector and its covariance
matrix of the time invariant Kalman filters are computed under the stationarity
condition:
```
z0  = (I - TTT)\CCC
P0 = reshape(I - kron(TTT, TTT))\vec(V), Nz, Nz)
```

Where:

- `kron(TTT, TTT)` is a matrix of dimension `Nz^2` x `Nz^2`, the Kronecker
  product of `TTT`
- `vec(V)` is the `Nz^2` x 1 column vector constructed by stacking the `Nz`
  columns of `V`

All eigenvalues of `TTT` are inside the unit circle when the state space model
is stationary.  When the preceding formula cannot be applied, the initial state
vector estimate is set to `CCC` and its covariance matrix is given by `1e6 * I`.
"""
function kalman_filter{S<:AbstractFloat}(data::Matrix{S},
    TTT::Matrix{S}, RRR::Matrix{S}, CCC::Vector{S},
    QQ::Matrix{S}, ZZ::Matrix{S}, DD::Vector{S}, MM::Matrix{S}, EE::Matrix{S},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = false, n_presample_periods::Int = 0)

    T = size(data, 2)
    regime_indices = Range{Int64}[1:T]

    kalman_filter(regime_indices, data, Matrix{S}[TTT], Matrix{S}[RRR], Vector{S}[CCC],
        Matrix{S}[QQ], Matrix{S}[ZZ], Vector{S}[DD],
        Matrix{S}[MM], Matrix{S}[EE], z0, P0;
        allout = allout, n_presample_periods = n_presample_periods)
end

function kalman_filter{S<:AbstractFloat}(regime_indices::Vector{Range{Int64}},
    data::Matrix{S}, TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}},
    QQs::Vector{Matrix{S}}, ZZs::Vector{Matrix{S}}, DDs::Vector{Vector{S}},
    MMs::Vector{Matrix{S}}, EEs::Vector{Matrix{S}},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = false, n_presample_periods::Int = 0)

    n_regimes = length(regime_indices)

    # Dimensions
    T  = size(data,    2) # number of periods of data
    Nz = size(TTTs[1], 1) # number of states
    Ne = size(RRRs[1], 2) # number of shocks
    Ny = size(ZZs[1],  1) # number of observables

    # Populate initial conditions if they are empty
    if isempty(z0) || isempty(P0)
        e, _ = eig(TTTs[1])
        if all(abs(e) .< 1.)
            z0 = (UniformScaling(1) - TTTs[1])\CCCs[1]
            P0 = solve_discrete_lyapunov(TTTs[1], RRRs[1]*QQs[1]*RRRs[1]')
        else
            z0 = CCC
            P0 = 1e6 * eye(Nz)
        end
    end

    z = z0
    P = P0

    # Initialize outputs
    if allout
        pred          = zeros(S, Nz, T)
        vpred         = zeros(S, Nz, Nz, T)
        yprederror    = NaN*zeros(S, Ny, T)
        ystdprederror = NaN*zeros(S, Ny, T)
        filt          = zeros(Nz, T)
        vfilt         = zeros(Nz, Nz, T)
    end

    log_likelihood = 0.0

    for i = 1:n_regimes
        # Get state-space system matrices for this regime
        regime_periods = regime_indices[i]
        regime_data = data[:, regime_periods]

        TTT, RRR, CCC = TTTs[i], RRRs[i], CCCs[i]
        QQ,  ZZ,  DD  = QQs[i],  ZZs[i],  DDs[i]
        MM,  EE       = MMs[i],  EEs[i]

        V = RRR*QQ*RRR'    # V = Var(z_t) = Var(Rϵ_t)
        R = EE + MM*QQ*MM' # R = Var(y_t) = Var(u_t)
        G = RRR*QQ*MM'     # G = Cov(z_t, y_t)

        for t in regime_periods
            # If an element of the vector y_t is missing (NaN) for the observation t, the
            # corresponding row is ditched from the measurement equation
            nonmissing = !isnan(data[:, t])
            y_t  = data[nonmissing, t]
            ZZ_t = ZZ[nonmissing, :]
            G_t  = G[:, nonmissing]
            R_t  = R[nonmissing, nonmissing]
            Ny_t = length(y_t)
            DD_t = DD[nonmissing]

            ## Forecast
            z = CCC + TTT*z                    # z_{t|t-1} = CCC + TTT*z_{t-1|t-1}
            P = TTT*P*TTT' + V                 # P_{t|t-1} = TTT*P_{t-1|t-1}*TTT' + TTT*Var(η_t)*TTT'
            dy = y_t - ZZ_t*z - DD_t           # dy = y_t - ZZ*z_{t|t-1} - DD is prediction error or innovation
            ZG = ZZ_t*G_t                      # ZG is ZZ*Cov(η_t, ϵ_t)
            D = ZZ_t*P*ZZ_t' + ZG + ZG' + R_t  # D = ZZ*P_{t|t-1}*ZZ' + ZG + ZG' + R_t
            D = (D+D')/2

            if allout
                pred[:, t]                   = z
                vpred[:, :, t]               = P
                yprederror[nonmissing, t]    = dy
                ystdprederror[nonmissing, t] = dy ./ sqrt(diag(D))
            end

            ddy = D\dy

            # We evaluate the log likelihood function by adding values of L at every iteration
            # step (for each t = 1,2,...T)
            if t > n_presample_periods
                log_likelihood += -log(det(D))/2 - first(dy'*ddy/2) - Ny_t*log(2*pi)/2
            end

            ## Update
            PZG = P*ZZ_t' + G_t
            z = z + PZG*ddy                    # z_{t|t} = z_{t|t-1} + P_{t|t-1}*ZZ(Θ)' + ...
            P = P - PZG/D*PZG'                 # P_{t|t} = P_{t|t-1} - PZG*(1/D)*PZG

            if allout
                filt[:, t]     = z
                vfilt[:, :, t] = P
            end

        end # of loop through this regime's periods

    end # of loop through regimes

    zend = z
    Pend = P

    if allout && n_presample_periods > 0
        mainsample_periods = n_presample_periods+1:T

        # If we choose to discard presample periods, then we reassign `z0`
        # and `P0` to be their values at the end of the presample/beginning
        # of the main sample
        z0 = squeeze(filt[:,     n_presample_periods], 2)
        P0 = squeeze(vfilt[:, :, n_presample_periods], 3)

        pred          = pred[:,     mainsample_periods]
        vpred         = vpred[:, :, mainsample_periods]
        yprederror    = yprederror[:,       mainsample_periods]
        ystdprederror = ypredstderror[:, :, mainsample_periods]
        filt          = filt[:,      mainsample_periods]
        vfilt         = vfilt[:, :, mainsample_periods]
    end

    if allout
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))

        return log_likelihood, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt, z0, P0
    else
        return log_likelihood
    end
end

"""
```
kalman_filter{S<:AbstractFloat}(m, data,
    TTT = Matrix{S}(0, 0), RRR = Matrix{S}(0, 0), CCC = Vector{S}(0,),
    z0 = Vector{S}(0,), P0 = Matrix{S}(0, 0);
    ZZ = Matrix{S}(0, 0), DD = Vector{S}(0,), QQ = Matrix{S}(0, 0),
    MM = Matrix{S}(0, 0), EE = Matrix{S}(0, 0), VVall = Matrix{S}(0, 0),
    allout = false, catch_errors = false,
    include_presample = true)
```

Implements the Kalman filter, accounting for the zero lower bound.

### Inputs

- `m`: model object
- `data`: `Ny` x `T` matrix containing data `y(1), ... , y(T)`.
- `TTT`: optional `Nz` x `Nz` matrix for a time-invariant transition matrix in
  the transition equation. If not provided, it will be calculated from `m`
- `RRR`: optional `Nz` x `Ne` matrix mapping exogenous shocks to states in the
  transition equation. If not provided, it will be calculated from `m`
- `CCC`: `Nz` x 1 vector, the constant term in the transition equation. If
  not provided, it will be calculated from `m`
- `z0`: optional `Nz` x 1 initial state vector
- `P0`: optional `Nz` x `Nz` covariance matrix of the initial state vector

where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `T`: number of periods of data

#### Keyword arguments

- `ZZ`: `Ny` x `Nz` matrix mapping states to observables
- `DD`: `Ny` x 1 vector, the constant term in the measurement equation
- `QQ`: `Ne` x `Ne` matrix of exogenous shock covariances
- `MM`: `Ny` x `Ne` matrix. See `Measurement` type for description
- `EE`: `Ny` x `Ny` matrix. See `Measurement` type for description
- `allout`: indicates whether we want optional output variables returned as well
- `include_presample`: indicates whether to include presample periods in the
  returned Kalman object. If true, we concatenate Kalman objects from all three
  regimes, else only R2 and R3

### Outputs

- a `Kalman` object. See documentation for `Kalman`
"""
function kalman_filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    catch_errors::Bool = false, allout::Bool = false, include_presample::Bool = true)

    # Compute the transition equation:
    #   S_t = CCC + TTT*S_{t-1} + RRR*ε_t
    # where
    #   Var(ε_t) = QQ
    # If we are in Metropolis-Hastings, then any errors coming out of `gensys`
    # should be caught and a -Inf posterior should be returned.
    TTT, RRR, CCC = Matrix{S}(), Matrix{S}(), Vector{S}()

    try
        TTT, RRR, CCC = solve(m)
    catch err
        if catch_errors && isa(err, GensysError)
            info(err.msg)
            return Kalman(-Inf, Vector{S}(), Matrix{S}())
        else
            rethrow(err)
        end
    end

    # Define the measurement equation:
    #   Y_t = DD + ZZ*S_t + u_t
    # where
    #   u_t = η_t + MM*ε_t
    #   Var(η_t) = EE
    meas = measurement(m, TTT, RRR, CCC; shocks = true)
    QQ, ZZ, DD = meas[:QQ], meas[:ZZ], meas[:DD]
    MM, EE     = meas[:MM], meas[:EE]

    # Compute the likelihood using the Kalman filter
    kalman_filter(m, data, TTT, RRR, CCC, QQ, ZZ, DD, MM, EE, z0, P0;
        allout = allout, include_presample = include_presample)
end

function kalman_filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    TTT::Matrix{S}, RRR::Matrix{S}, CCC::Vector{S},
    QQ::Matrix{S}, ZZ::Matrix{S}, DD::Vector{S}, MM::Matrix{S}, EE::Matrix{S},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = false, include_presample::Bool = true)

    # If z0 and P0 provided, check that rows and columns corresponding to
    # anticipated shocks are zero in P0
    if !isempty(z0) && !isempty(P0)
        ant_state_inds = setdiff(1:n_states_augmented(m), inds_states_no_ant(m))
        @assert all(x -> x == 0, P0[:, ant_state_inds])
        @assert all(x -> x == 0, P0[ant_state_inds, :])
    end

    # Partition sample into pre- and post-ZLB regimes
    # Note that the post-ZLB regime may be empty if we do not impose the ZLB
    T = size(data, 2)

    regime_indices = Vector{Range{Int64}}(2)
    if n_anticipated_shocks(m) > 0
        regime_indices[1] = 1:index_zlb_start(m)-1
        regime_indices[2] = index_zlb_start(m):T # allows for conditional data
    else
        regime_indices[1] = 1:T
        regime_indices[2] = 1:0
    end

    # For the pre-ZLB periods, we must zero out the shocks corresponding to
    # anticipated shocks
    QQ_preZLB = copy(QQ)

    ant_shock_inds = setdiff(1:n_shocks_exogenous(m), inds_shocks_no_ant(m))
    for i in ant_shock_inds
        QQ_preZLB[i, i] = 0
    end

    # Specify number of presample periods if we don't want to include them in
    # the final results
    T0 = include_presample ? 0 : n_presample_periods(m)

    # Run Kalman filter, construct Kalman object, and return
    out = kalman_filter(regime_indices, data, fill(TTT, 2), fill(RRR, 2), fill(CCC, 2),
              Matrix{S}[QQ_preZLB, QQ], fill(ZZ, 2), fill(DD, 2), fill(MM, 2), fill(EE, 2),
              z0, P0; allout = allout, n_presample_periods = T0)

    return Kalman(out...)
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
- `P0`: variance-covariance matrix for `z0`

#### Fields filled in when `allout=true` in a call to `kalman_filter`:
- `pred`: a `Nz` x `T` matrix containing one-step predicted state vectors
- `vpred`: a `Nz` x `Nz` x `T` matrix containing mean square errors of predicted
  state vectors

- `filt`: an `Nz` x `T` matrix containing filtered state vectors
- `vfilt`: an `Nz` x `Nz` x `T` matrix containing mean square errors of filtered
  state vectors
"""
immutable Kalman{S<:AbstractFloat}
    L::S                  # likelihood
    zend::Vector{S}       # last-period state vector
    Pend::Matrix{S}       # last-period variance-covariance matrix for the states
    pred::Matrix{S}       # predicted value of states in period T+1
    vpred::Array{S, 3}    # predicted variance-covariance matrix for states in period T+1
    yprederror::Matrix{S}
    ystdprederror::Matrix{S}
    rmse::Matrix{S}
    rmsd::Matrix{S}
    filt::Matrix{S}        # filtered states
    vfilt::Array{S, 3}     # mean squared errors of filtered state vectors
    z0::Vector{S}          # starting-period state vector
    vz0::Matrix{S}         # starting-period variance-covariance matrix for the states
end

function Kalman{S<:AbstractFloat}(L::S,
                                  zend::Vector{S}          = Vector{S}(),
                                  Pend::Matrix{S}          = Matrix{S}(),
                                  pred::Matrix{S}          = Matrix{S}(),
                                  vpred::Array{S, 3}       = Array{S}(0, 0, 0);
                                  yprederror::Matrix{S}    = Matrix{S}(),
                                  ystdprederror::Matrix{S} = Matrix{S}(),
                                  rmse::Matrix{S}          = Matrix{S}(),
                                  rmsd::Matrix{S}          = Matrix{S}(),
                                  filt::Matrix{S}          = Matrix{S}(),
                                  vfilt::Array{S, 3}       = Array{S}(0, 0, 0),
                                  z0::Vector{S}            = Vector{S}(),
                                  P0::Matrix{S}            = Matrix{S}())
    return Kalman{S}(L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt, z0, P0)
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
    k2::Kalman{S}; allout::Bool = false)

    L = k1[:L] + k2[:L]
    zend = k2[:zend]
    Pend = k2[:Pend]

    if allout
        pred  = hcat(k1[:pred], k2[:pred])
        vpred = cat(3, k1[:vpred], k2[:vpred])
        yprederror    = hcat(k1[:yprederror], k2[:yprederror])
        ystdprederror = hcat(k1[:ystdprederror], k2[:yprederror])
        rmse  = sqrt(mean((yprederror.^2)', 1))
        rmsd  = sqrt(mean((ystdprederror.^2)', 1))
        filt  = hcat(k1[:filt], k2[:filt])
        vfilt = cat(3, k1[:vfilt], k2[:vfilt])
        z0    = k1[:z0]
        P0    = k1[:vz0]

        return Kalman(L, zend, Pend, pred, vpred, yprederror, ystdprederror,
            rmse, rmsd, filt, vfilt, z0, P0)
    else
        return Kalman(L, zend, Pend)
    end
end