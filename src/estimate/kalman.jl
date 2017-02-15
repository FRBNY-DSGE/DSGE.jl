#=
This code is loosely based on a routine originally copyright Federal Reserve Bank of Atlanta
and written by Iskander Karibzhanov.
=#

"""
```
kalman_filter{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    TTT::Matrix{S}, CCC::Vector{S}, ZZ::Matrix{S}, DD::Vector{S}, VVall::Matrix{S},
    z0::Vector{S} = Vector{S}(), vz0::Matrix{S} = Matrix{S}(); lead::Int = 0,
    allout::Bool = false, include_presample::Bool = true)
```

### Inputs

- `m`: model object
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
- `vz0`: an optional `Nz` x `Nz` covariance matrix of an initial state vector.
- `lead`: the number of steps to forecast after the end of the data.
- `allout`: an optional keyword argument indicating whether we want optional output
  variables returned as well
- `include_presample`: indicates whether to include presample periods in the
  returned Kalman object. If `!include_presample`, then we don't add presample
  periods to the likelohood, and we also set the `z0` and `vz0` fields in the
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

When `z0` and `Vz0` are omitted, the initial state vector and its covariance
matrix of the time invariant Kalman filters are computed under the stationarity
condition:
```
z0  = (I - TTT)\CCC
vz0 = reshape(I - kron(TTT, TTT))\vec(V), Nz, Nz)
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
function kalman_filter{S<:AbstractFloat}(m::AbstractModel,
                                         data::Matrix{S},
                                         TTT::Matrix{S},
                                         CCC::Vector{S},
                                         ZZ::Matrix{S},
                                         DD::Vector{S},
                                         VVall::Matrix{S},
                                         z0::Vector{S} = Vector{S}(),
                                         vz0::Matrix{S} = Matrix{S}();
                                         lead::Int = 0,
                                         allout::Bool = false,
                                         include_presample::Bool = true)
    T = size(data, 2)
    Nz = length(CCC)
    Ny = length(DD)
    V = VVall[1:Nz, 1:Nz] # V = RRR*QQ*RRR'

    if isempty(z0) || isempty(vz0)
        e, _ = eig(TTT)
        if all(abs(e) .< 1.)
            z0  = (eye(Nz) - TTT)\CCC
            vz0 = solve_discrete_lyapunov(TTT, V)
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
        DD_t = DD[nonmissing]              # DD_t


        ## forecasting
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
        #   step (for each t = 1,2,...T)
        if include_presample || (!include_presample && t > n_presample_periods(m))
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

        # If !include_presample, then we reassign `z0` and `P0` to be their
        # values at the end of the presample/beginning of the main sample
        if !include_presample && t == n_presample_periods(m)
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
kalman_filter_2part{S<:AbstractFloat}(m, data,
    TTT = Matrix{S}(0, 0), RRR = Matrix{S}(0, 0), CCC = Vector{S}(0,),
    z0 = Vector{S}(0,), vz0 = Matrix{S}(0, 0);
    ZZ = Matrix{S}(0, 0), DD = Vector{S}(0,), QQ = Matrix{S}(0, 0),
    MM = Matrix{S}(0, 0), EE = Matrix{S}(0, 0), VVall = Matrix{S}(0, 0),
    lead = 0, allout = false, catch_errors = false,
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
- `vz0`: optional `Nz` x `Nz` covariance matrix of the initial state vector

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
- `VVall`: `Nz + Ny` x `Nz + Ny` matrix. See `Measurement` type for description
- `lead`: number of steps to forecast after the end of the data
- `allout`: indicates whether we want optional output variables returned as well
- `include_presample`: indicates whether to include presample periods in the
  returned Kalman object. If true, we concatenate Kalman objects from all three
  regimes, else only R2 and R3

### Outputs

- a `Kalman` object. See documentation for `Kalman`
"""
function kalman_filter_2part{S<:AbstractFloat}(m::AbstractModel,
                                               data::Matrix{S},
                                               TTT::Matrix{S}   = Matrix{S}(0, 0),
                                               RRR::Matrix{S}   = Matrix{S}(0, 0),
                                               CCC::Vector{S}   = Vector{S}(0,),
                                               z0::Vector{S}    = Vector{S}(0,),
                                               vz0::Matrix{S}   = Matrix{S}(0, 0);
                                               ZZ::Matrix{S}    = Matrix{S}(0, 0),
                                               DD::Vector{S}    = Vector{S}(0,),
                                               QQ::Matrix{S}    = Matrix{S}(0, 0),
                                               MM::Matrix{S}    = Matrix{S}(0, 0),
                                               EE::Matrix{S}    = Matrix{S}(0, 0),
                                               VVall::Matrix{S} = Matrix{S}(0, 0),
                                               lead::Int        = 0,
                                               allout::Bool     = false,
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
                return LIKE_NULL_OUTPUT
            else
            rethrow(err)
            end
        end
    end

    # Step 2: Define the measurement equation:
    #   Y_t = DD + ZZ*S_t + u_t
    # where
    #   Var(η_t) = EE
    #   u_t = η_t + MM*ε_t
    #   Var(u_t) = HH = EE+MM QQ MM'
    #   Cov(ε_t,u_t) = VV = QQ*MM'
    if isempty(ZZ) || isempty(DD) || isempty(QQ) || isempty(MM) || isempty(EE) || isempty(VVall)
        meas = measurement(m, TTT, RRR, CCC; shocks = true)
        ZZ, QQ, VVall = meas[:ZZ], meas[:QQ], meas[:VVall]
        MM, EE = meas[:MM], meas[:EE]

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
    R2[:VVall] = [[RRR*R2[:QQ]*RRR'  RRR*R2[:QQ]*MM'];
                  [MM*R2[:QQ]*RRR'   EE+MM*R2[:QQ]*MM']]

    for d in [:QQ, :VVall]
        R1[d] = R2[d]
    end

    # Step 3: Compute log-likelihood using the Kalman filter.
    # Note that `kalman_filter` assumes a transition equation of the form:
    #   S_t = TTT*S_{t-1} + ε2_t
    # where ε2_t = RRR*ε_t. Therefore redefine:
    #   QQ2 = Var(ε2_t) = RRR*QQ*RRR'
    #   VV2 = Cov(ε2_t, u_t) = RRR*VV
    #   VVall = Var([ε2_t; u_t])    (joint variance of the two shocks)

    # Run Kalman filter on presample, calculating `z0` and `vz0` in
    # `kalman_filter` if necessary
    if isempty(z0) || isempty(vz0)
        k1 = kalman_filter(m, R1[:data], TTT, CCC, ZZ, DD, R1[:VVall];
            allout = allout, include_presample = true)
    else
        # Check that rows and columns corresponding to anticipated shocks are zero in vz0
        ant_state_inds = setdiff(1:n_states_augmented(m), inds_states_no_ant(m))
        @assert all(x -> x == 0, vz0[:, ant_state_inds])
        @assert all(x -> x == 0, vz0[ant_state_inds, :])

        k1 = kalman_filter(m, R1[:data], TTT, CCC, ZZ, DD, R1[:VVall], z0, vz0;
                           allout = allout, include_presample = true)
    end

    # Run Kalman filter on normal period
    k2 = kalman_filter(m, R2[:data], TTT, CCC, ZZ, DD, R2[:VVall], k1[:zend], k1[:Pend];
                       allout = allout, include_presample = true)

    # Run Kalman filter on ZLB period
    k3 = kalman_filter(m, R3[:data], TTT, CCC, ZZ, DD, VVall, k2[:zend], k2[:Pend];
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
        vz0   = k1[:vz0]

        return Kalman(L, zend, Pend, pred, vpred, yprederror, ystdprederror,
            rmse, rmsd, filt, vfilt, z0, vz0)
    else
        return Kalman(L, zend, Pend)
    end
end