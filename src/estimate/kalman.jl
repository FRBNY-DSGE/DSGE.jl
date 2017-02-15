#=
This code is loosely based on a routine originally copyright Federal Reserve Bank of Atlanta
and written by Iskander Karibzhanov.
=#

"""
```
kalman_filter{S<:AbstractFloat}(data::Matrix{S},
    TTT::Matrix{S}, CCC::Vector{S}, ZZ::Matrix{S}, DD::Vector{S}, VVall::Matrix{S},
    z0::Vector{S} = Vector{S}(), P0::Matrix{S} = Matrix{S}();
    allout::Bool = false, include_presample::Bool = true)
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
- `include_presample`: indicates whether to include presample periods in the
  returned Kalman object. If `!include_presample`, then we don't add presample
  periods to the likelohood, and we also set the `z0` and `P0` fields in the
  returned `Kalman` object to be the states and variance-covariance matrices at
  the end of the presample/beginning of the main sample

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
    allout::Bool = false, include_presample::Bool = true)

    # Dimensions
    T  = size(data, 2) # number of periods of data
    Nz = size(TTT, 1)  # number of states
    Ne = size(RRR, 2)  # number of shocks
    Ny = size(ZZ, 1)   # number of observables

    # Populate initial conditions if they are empty
    if isempty(z0) || isempty(P0)
        e, _ = eig(TTT)
        if all(abs(e) .< 1.)
            z0 = (UniformScaling(1) - TTT)\CCC
            P0 = solve_discrete_lyapunov(TTT, RRR*QQ*RRR')
        else
            z0 = CCC
            P0 = 1e6 * eye(Nz)
        end
    end

    z = z0
    P = P0

    # Check input matrix dimensions
    @assert size(data, 1) == Ny
    @assert size(TTT) == (Nz, Nz)
    @assert size(RRR) == (Nz, Ne)
    @assert size(CCC) == (Nz,)
    @assert size(QQ)  == (Ne, Ne)
    @assert size(ZZ)  == (Ny, Nz)
    @assert size(DD)  == (Ny,)
    @assert size(MM)  == (Ny, Ne)
    @assert size(EE)  == (Ny, Ny)
    @assert size(z)   == (Nz,)
    @assert size(P)   == (Nz, Nz)

    if allout
        pred          = zeros(S, Nz, T)
        vpred         = zeros(S, Nz, Nz, T)
        yprederror    = NaN*zeros(S, Ny, T)
        ystdprederror = NaN*zeros(S, Ny, T)
        filt          = zeros(Nz, T)
        vfilt         = zeros(Nz, Nz, T)
    end

    V = RRR*QQ*RRR'    # V = Var(z_t) = Var(Rϵ_t)
    R = EE + MM*QQ*MM' # R = Var(y_t) = Var(u_t)
    G = RRR*QQ*MM'     # G = Cov(z_t, y_t)

    L = zero(S)

    for t = 1:T
        # If an element of the vector y(t) is missing (NaN) for the observation t, the
        # corresponding row is ditched from the measurement equation
        nonmissing = !isnan(data[:, t])
        data_t = data[nonmissing, t]
        ZZ_t = ZZ[nonmissing, :]
        G_t  = G[:, nonmissing]
        R_t  = R[nonmissing, nonmissing]
        Ny_t = length(data_t)
        DD_t = DD[nonmissing]

        ## Forecast
        z = CCC + TTT*z                    # z_{t|t-1} = CCC + TTT(Θ)*z_{t-1|t-1}
        P = TTT*P*TTT' + V                 # P_{t|t-1} = TTT(Θ)*P_{t-1|t-1}*TTT(Θ)' + TTT(Θ)*Var(η_t)*TTT(Θ)'
        dy = data_t - ZZ_t*z - DD_t        # dy = y_t - ZZ(Θ)*z_{t|t-1} - DD is prediction error or innovation
        ZG = ZZ_t*G_t                      # ZG is ZZ*Cov(η_t, ϵ_t)
        D = ZZ_t*P*ZZ_t' + ZG + ZG' + R_t  # D = ZZ*P_{t|t-1}*ZZ' + ZG + ZG' + R_t
        D = (D+D')/2

        if allout
            pred[:, t]                   = z
            vpred[:, :, t]               = P
            yprederror[nonmissing, t]    = dy
            ystdprederror[nonmissing, t] = dy./sqrt(diag(D))
        end

        ddy = D\dy

        # We evaluate the log likelihood function by adding values of L at every iteration
        # step (for each t = 1,2,...T)
        if include_presample || (!include_presample && t > n_presample_periods(m))
            L += -log(det(D))/2 - first(dy'*ddy/2) - Ny_t*log(2*pi)/2
        end

        ## Update
        PZG = P*ZZ_t' + G_t
        z = z + PZG*ddy                    # z_{t|t} = z_{t|t-1} + P_{t|t-1}*ZZ(Θ)' + ...
        P = P - PZG/D*PZG'                 # P_{t|t} = P_{t|t-1} - PZG*(1/D)*PZG

        if allout
            PZZ = P*ZZ_t'
            filt[:, t]     = z
            vfilt[:, :, t] = P
        end

        # If !include_presample, then we reassign `z0` and `P0` to be their
        # values at the end of the presample/beginning of the main sample
        if !include_presample && t == n_presample_periods(m)
            z0 = z
            P0 = P
        end
    end

    zend = z
    Pend = P

    if allout
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))

        return Kalman(L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt, z0, P0)
    else
        return Kalman(L, zend, Pend)
    end

end

"""
```
kalman_filter_2part{S<:AbstractFloat}(m, data,
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
function kalman_filter_2part{S<:AbstractFloat}(m::AbstractModel,
                                               data::Matrix{S},
                                               TTT::Matrix{S} = Matrix{S}(0, 0),
                                               RRR::Matrix{S} = Matrix{S}(0, 0),
                                               CCC::Vector{S} = Vector{S}(0,),
                                               z0::Vector{S}  = Vector{S}(0,),
                                               P0::Matrix{S}  = Matrix{S}(0, 0);
                                               QQ::Matrix{S}  = Matrix{S}(0, 0),
                                               ZZ::Matrix{S}  = Matrix{S}(0, 0),
                                               DD::Vector{S}  = Vector{S}(0,),
                                               MM::Matrix{S}  = Matrix{S}(0, 0),
                                               EE::Matrix{S}  = Matrix{S}(0, 0),
                                               allout::Bool   = false,
                                               catch_errors::Bool = false,
                                               include_presample::Bool = true)

    # Partition sample into three regimes, and store associated matrices:
    # - R1: presample
    # - R2: normal
    # - R3: zero lower bound and beyond
    R1 = Dict{Symbol, Array{S}}()
    R2 = Dict{Symbol, Array{S}}()
    R3 = Dict{Symbol, Array{S}}()

    R1[:data] = data[:, inds_presample_periods(m)]
    R2[:data] = data[:, inds_prezlb_periods(m)]
    R3[:data] = data[:, index_zlb_start(m):end] # allows for conditional data

    # Step 1: Compute the transition equation:
    #   S_t = CCC + TTT*S_{t-1} + RRR*ε_t
    # where
    #   Var(ε_t) = QQ
    # If we are in Metropolis-Hastings, then any errors coming out of `gensys`
    # should be caught and a -Inf posterior should be returned.
    if isempty(TTT) || isempty(RRR) || isempty(CCC)
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
    end

    # Step 2: Define the measurement equation:
    #   Y_t = DD + ZZ*S_t + u_t
    # where
    #   u_t = η_t + MM*ε_t
    #   Var(η_t) = EE
    if isempty(QQ) || isempty(ZZ) || isempty(DD) || isempty(MM) || isempty(EE)
        meas = measurement(m, TTT, RRR, CCC; shocks = true)
        QQ, ZZ, MM, EE = meas[:QQ], meas[:ZZ], meas[:MM], meas[:EE]

        # If DD specifically is nonempty (most often a vector of zeros, as in
        # the Durbin-Koopman smoother), we want to use that DD instead of the
        # one calculated from the measurement equation
        DD = isempty(DD) ? meas[:DD] : DD
    end

    # For the pre-ZLB periods, we must zero out the shocks corresponding to anticipated shocks
    R2[:QQ] = copy(QQ)
    ant_shock_inds = setdiff(1:n_shocks_exogenous(m), inds_shocks_no_ant(m))
    for i in ant_shock_inds
        R2[:QQ][i, i] = 0
    end

    R1[:QQ] = R2[:QQ]

    # Step 3: Compute log-likelihood using the Kalman filter

    # Run Kalman filter on presample, calculating `z0` and `P0` in
    # `kalman_filter` if necessary
    if isempty(z0) || isempty(P0)
        k1 = kalman_filter(m, R1[:data], TTT, RRR, CCC, R1[:QQ], ZZ, DD, MM, EE,
            allout = allout, include_presample = true)
    else
        # Check that rows and columns corresponding to anticipated shocks are zero in P0
        ant_state_inds = setdiff(1:n_states_augmented(m), inds_states_no_ant(m))
        @assert all(x -> x == 0, P0[:, ant_state_inds])
        @assert all(x -> x == 0, P0[ant_state_inds, :])

        k1 = kalman_filter(m, R1[:data], TTT, RRR, CCC, R1[:QQ], ZZ, DD, MM, EE, z0, P0;
                           allout = allout, include_presample = true)
    end

    # Run Kalman filter on normal period
    k2 = kalman_filter(m, R2[:data], TTT, RRR, CCC, R2[:QQ], ZZ, DD, MM, EE, k1[:zend], k1[:Pend];
                       allout = allout, include_presample = true)

    # Run Kalman filter on ZLB period
    k3 = kalman_filter(m, R3[:data], TTT, RRR, CCC, QQ, ZZ, DD, MM, EE, k2[:zend], k2[:Pend];
                       allout = allout, include_presample = true)

    # Concatenate Kalman objects
    if include_presample
        k12 = cat(m, k1, k2; allout = allout)
        kal = cat(m, k12, k3; allout = allout)
    else
        kal = cat(m, k2, k3; allout = allout)
    end

    ## Return concatenated Kalman
    return kal
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
                                  P0::Matrix{S}           = Matrix{S}(0, 0))
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