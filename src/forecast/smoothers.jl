"""
```
kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    cond_type::Symbol = :none, include_presample::Bool = false)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    include_presample::Bool = false)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    cond_type::Symbol = :none, include_presample::Bool = false)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    include_presample::Bool = false)
```
This is a Kalman Smoothing program based on S.J. Koopman's \"Disturbance
Smoother for State Space Models\" (Biometrika, 1993), as specified in
Durbin and Koopman's \"A Simple and Efficient Simulation Smoother for
State Space Time Series Analysis\" (Biometrika, 2002). The algorithm has been
simplified for the case in which there is no measurement error, and the
model matrices do not vary with time.

Unlike other Kalman Smoothing programs, there is no need to invert
singular matrices using the Moore-Penrose pseudoinverse (pinv), which
should lead to efficiency gains and fewer inversion problems. Also, the
states vector and the corresponding matrices do not need to be augmented
to include the shock innovations. Instead they are saved automatically
in the `eta_hat` matrix.

### Inputs:

- `m`: model object
- `data`: the (`Ny` x `Nt`) matrix of observable data
- `T`: the (`Nz` x `Nz`) transition matrix
- `R`: the (`Nz` x `Ne`) matrix translating shocks to states
- `C`: the (`Nz` x 1) constant vector in the transition equation
- `Q`: the (`Ne` x `Ne`) covariance matrix for the shocks
- `Z`: the (`Ny` x `Nz`) measurement matrix
- `D`: the (`Ny` x 1) constant vector in the measurement equation
- `A0`: the (`Nz` x 1) initial (time 0) states vector
- `P0`: the (`Nz` x `Nz`) initial (time 0) state covariance matrix
- `pred`: the (`Nz` x `Nt`) matrix of one-step-ahead predicted states (from the
  Kalman Filter)
- `vpred`: the (`Nz` x `Nz` x `Nt`) matrix of one-step-ahead predicted
  covariance matrices
- `cond_type`: optional keyword argument specifying the conditional data type:
  one of `:none`, `:semi`, or `:full`. This is only necessary when a DataFrame
  (as opposed to a data matrix) is passed in, so that `df_to_matrix` knows how
  many periods of data to keep
- `include_presample`: indicates whether or not to return presample periods in
  the returned smoothed states and shocks. Defaults to `false`

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

- `α_hat`: the (`Nz` x `Nt`) matrix of smoothed states
- `η_hat`: the (`Ne` x `Nt`) matrix of smoothed shocks

If `n_presample_periods(m)` is nonzero, the `α_hat` and `η_hat` matrices will be
shorter by that number of columns (taken from the beginning).

### Notes

The state space model is defined as follows:
```
y(t) = Z*α(t) + D             (state or transition equation)
α(t+1) = T*α(t) + R*η(t+1)    (measurement or observation equation)
```
"""
function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # extract system matrices
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]

    # call actual Kalman smoother
    kalman_smoother(m, df, T, R, C, Q, Z, D, A0, P0, pred, vpred; cond_type =
        cond_type, include_presample = include_presample)
end

function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    include_presample::Bool = false)

    # extract system matrices
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]

    # call actual Kalman smoother
    kalman_smoother(m, data, T, R, C, Q, Z, D, A0, P0, pred, vpred;
        include_presample = include_presample)
end

function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # convert DataFrame to matrix
    data = df_to_matrix(m, df; cond_type = cond_type)

    # call actual Kalman smoother
    kalman_smoother(m, data, T, R, C, Q, Z, D, A0, P0, pred, vpred;
        include_presample = include_presample)
end

function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    include_presample::Bool = false)

    Ne = size(R, 2)
    Ny = size(data, 1)
    Nt = size(data, 2)
    Nz = size(T, 1)

    # Check data is well-formed wrt model settings
    @assert Ny == n_observables(m)
    @assert Nt >= n_presample_periods(m) + n_prezlb_periods(m) + n_zlb_periods(m)

    # Anticipated monetary policy shocks
    n_ant_shocks = n_anticipated_shocks(m)
    t_zlb_start  = index_zlb_start(m)

    r, η_hat = disturbance_smoother(m, data, T, R, C, Q, Z, D, pred, vpred)

    α_hat = zeros(Nz, Nt)
    ah_t = A0 + P0*r[:, 1]
    α_hat[:, 1] = ah_t

    shock_inds = inds_shocks_no_ant(m)
    for t = 2:Nt

        # This section relates to the zero bound framework, in which no
        # anticipated shocks are supposed to occur before the model switch.  In
        # these periods, this is accomplished by setting the relevant rows and
        # columns of the Q matrix to zero. In other periods, or in
        # specifications with zero bound off (and hence with
        # n_anticipated_shocks(m) = 0), the normal Q matrix can be used.

        if n_ant_shocks > 0
            if t < t_zlb_start
                Q_t = zeros(Ne, Ne)
                Q_t[shock_inds, shock_inds] = Q[shock_inds, shock_inds]
                ah_t = T*ah_t + R*Q_t*R'*r[:, t]
            else
                ah_t = T*ah_t + R*Q*R'*r[:, t]
            end
        else
            ah_t = T*ah_t + R*Q*R'*r[:, t]
        end

        α_hat[:, t] = ah_t
    end

    if !include_presample
        α_hat = α_hat[:, index_mainsample_start(m):end]
        η_hat = η_hat[:, index_mainsample_start(m):end]
    end

    return α_hat, η_hat
end

"""
```
disturbance_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Array{S}, pred::Matrix{S}, vpred::Array{S, 3})
```

This is a Kalman Smoothing program based on S.J. Koopman's \"Disturbance
Smoother for State Space Models\" (Biometrika, 1993), as specified in
Durbin and Koopman's \"A Simple and Efficient Simulation Smoother for
State Space Time Series Analysis\" (Biometrika, 2002). The algorithm has been
simplified for the case in which there is no measurement error, and the
model matrices do not vary with time.

This disturbance smoother is intended for use with the state smoother
`kalman_smoother` from the same papers (Koopman 1993, Durbin and Koopman
2002). It produces a matrix of vectors, `r`, that is used for state
smoothing, and an optional matrix, `eta_hat`, containing the smoothed
shocks. It has been adjusted to account for the possibility of missing
values in the data, and to accommodate the zero bound model, which
requires that no anticipated shocks occur before the zero bound window,
which is achieved by setting the entries in the `Q` matrix corresponding to
the anticipated shocks to zero in those periods.

### Inputs:

- `m`: model object
- `data`: the (`Ny` x `Nt`) matrix of observable data
- `T`: the (`Nz` x `Nz`) transition matrix
- `R`: the (`Nz` x `Ne`) matrix translating shocks to states
- `C`: the (`Nz` x 1) constant vector in the transition equation
- `Q`: the (`Ne` x `Ne`) covariance matrix for the shocks
- `Z`: the (`Ny` x `Nz`) measurement matrix
- `D`: the (`Ny` x 1) constant vector in the measurement equation
- `pred`: the (`Nz` x `Nt`) matrix of one-step-ahead predicted states (from the
  Kalman Filter)
- `vpred`: the (`Nz` x `Nz` x `Nt`) matrix of one-step-ahead predicted
  covariance matrices

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

- `r`: the (`Nz` x `Nt`) matrix used for state smoothing
- `η_hat`: the (`Ne` x `Nt`) matrix of smoothed shocks

### Notes

The state space model is defined as follows:
```
y(t) = Z*α(t) + D             (state or transition equation)
α(t+1) = T*α(t) + R*η(t+1)    (measurement or observation equation)
```
"""
function disturbance_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Vector{S}, pred::Matrix{S}, vpred::Array{S, 3})

    Nt = size(data, 2)
    Nz = size(T, 1)

    r = zeros(Nz, Nt) # holds r_{T-1}, ..., r_0
    r_t = zeros(Nz, 1)

    Ne = size(R, 2)
    η_hat = zeros(Ne, Nt)

    # Anticipated policy shocks metadata
    n_ant_shocks = n_anticipated_shocks(m)
    t_zlb_start  = index_zlb_start(m)

    for t = Nt:-1:1
        data_t = data[:, t]

        # This section deals with the possibility of missing values in the y_t
        # vector (especially relevant for smoothing over peachdata).
        nonmissing = !isnan(data_t)
        data_t = data_t[nonmissing]
        Z_t = Z[nonmissing, :]
        D_t = D[nonmissing]

        a = pred[:, t]
        P = vpred[:, :, t]

        F = Z_t*P*Z_t'
        v = data_t - Z_t*a - D_t
        K = T*P*Z_t'/F
        L = T - K*Z_t

        r_t = Z_t'/F*v + L'*r_t
        r[:, t] = r_t

        # This section relates to the zero bound framework, in which no
        # anticipated shocks are supposed to occur before the model switch.  In
        # these periods, this is accomplished by setting the relevant rows and
        # columns of the Q matrix to zero. In other periods, or in
        # specifications with zero bound off (and hence with
        # n_anticipated_shocks(m) = 0), the normal Q matrix can be used.
        shock_inds = inds_shocks_no_ant(m)
        if n_ant_shocks > 0
            if t < t_zlb_start
                Q_t = zeros(Ne, Ne)
                Q_t[shock_inds, shock_inds] = Q[shock_inds, shock_inds]
                η_hat[:, t] = Q_t * R' * r_t
            else
                η_hat[:, t] = Q * R' * r_t
            end
        else
            η_hat[:, t] = Q * R' * r_t
        end
    end

    return r, η_hat
end

"""
```
durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, system::System, A0::Vector{S}, P0::Matrix{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data:Matrix{S}, system::System, A0::Vector{S}, P0::Matrix{S};
    include_presample::Bool = false)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, A0::Array{S}, P0::Matrix{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, A0::Array{S}, P0::Matrix{S};
    include_presample::Bool = false)
```
This program is a simulation smoother based on Durbin and Koopman's
\"A Simple and Efficient Simulation Smoother for State Space Time Series
Analysis\" (Biometrika, 2002). The algorithm has been simplified for the
case in which there is no measurement error, and the model matrices do
not vary with time.

Unlike other simulation smoothers (for example, that of Carter and Kohn,
1994), this method does not require separate draws for each period, draws
of the state vectors, or even draws from a conditional distribution.
Instead, vectors of shocks are drawn from the unconditional distribution
of shocks, which is then corrected (via a Kalman Smoothing step), to
yield a draw of shocks conditional on the data. This is then used to
generate a draw of states conditional on the data. Drawing the states in
this way is much more efficient than other methods, as it avoids the need
for multiple draws of state vectors (requiring singular value
decompositions), as well as inverting state covariance matrices
(requiring the use of the computationally intensive and relatively
erratic Moore-Penrose pseudoinverse).

### Inputs:

- `m`: model object
- `data`: the (`Ny` x `Nt`) matrix of observable data
- `T`: the (`Nz` x `Nz`) transition matrix
- `R`: the (`Nz` x `Ne`) matrix translating shocks to states
- `C`: the (`Nz` x 1) constant vector in the transition equation
- `Q`: the (`Ne` x `Ne`) covariance matrix for the shocks
- `Z`: the (`Ny` x `Nz`) measurement matrix
- `D`: the (`Ny` x 1) constant vector in the measurement equation
- `A0`: the (`Nz` x 1) initial (time 0) states vector
- `P0`: the (`Nz` x `Nz`) initial (time 0) state covariance matrix. If
- `cond_type`: optional keyword argument specifying the conditional data type:
  one of `:none`, `:semi`, or `:full`. This is only necessary when a DataFrame
  (as opposed to a data matrix) is passed in, so that `df_to_matrix` knows how
  many periods of data to keep
- `include_presample`: indicates whether or not to return presample periods in
  the returned smoothed states and shocks. Defaults to `false`

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

- `α_hat`: the (`Nz` x `Nt`) matrix of smoothed states.
- `η_hat`: the (`Ne` x `Nt`) matrix of smoothed shocks.

If `n_presample_periods(m)` is nonzero, the `α_hat` and `η_hat` matrices will be
shorter by that number of columns (taken from the beginning).

### Notes

The state space model is defined as follows:
```
y(t) = Z*α(t) + D             (state or transition equation)
α(t+1) = T*α(t) + R*η(t+1)    (measurement or observation equation)
```
"""
function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, system::System, A0::Vector{S}, P0::Matrix{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # extract system matrices
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]
    M, E, V_all = system[:MM], system[:EE], system[:VVall]

    # call actual Durbin-Koopman smoother
    durbin_koopman_smoother(m, df, T, R, C, Q, Z, D, M, E, V_all, A0, P0;
        cond_type = cond_type, include_presample = include_presample)
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, system::System, A0::Vector{S}, P0::Matrix{S};
    include_presample::Bool = false)

    # extract system matrices
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]
    M, E, V_all = system[:MM], system[:EE], system[:VVall]

    # call actual Durbin-Koopman smoother
    durbin_koopman_smoother(m, data, T, R, C, Q, Z, D, M, E, V_all, A0, P0;
        include_presample = include_presample)
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, T::Matrix{S}, R::Matrix{S}, C::Array{S},
    Q::Matrix{S}, Z::Matrix{S}, D::Vector{S},
    M::Matrix{S}, E::Matrix{S}, V_all::Matrix{S},
    A0::Vector{S}, P0::Matrix{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # convert DataFrame to Matrix
    data = df_to_matrix(m, df; cond_type = cond_type)

    # call actual simulation smoother
    durbin_koopman_smoother(m, data, T, R, C, Q, Z, D, M, E, V_all, A0, P0;
        include_presample = include_presample)
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S},
    Q::Matrix{S}, Z::Matrix{S}, D::Vector{S},
    M::Matrix{S}, E::Matrix{S}, V_all::Matrix{S},
    A0::Array{S}, P0::Matrix{S}; include_presample::Bool = false)

    # Get matrix dimensions
    Ny = size(data, 1)
    Nt = size(data, 2)
    Nz = size(T, 1)
    Ne = size(R, 2)

    # Check data is well-formed wrt model settings
    @assert Ny == n_observables(m)
    @assert Nt >= n_presample_periods(m) + n_prezlb_periods(m) + n_zlb_periods(m)

    # Anticipated monetary policy shocks
    n_ant_shocks = n_anticipated_shocks(m)
    t_zlb_start  = index_zlb_start(m)

    # Draw initial state α_0+ and sequence of shocks η+
    U, eig, _ = svd(P0)
    dist_α = DegenerateMvNormal(zeros(S, Nz), U * diagm(sqrt(eig)))
    dist_η = DegenerateMvNormal(zeros(S, Ne), sqrt(Q))

    if m.testing
        α_plus_0 = zeros(S, Nz)
        η_plus   = zeros(S, Ne, Nt)
    else
        α_plus_0 = rand(dist_α)
        η_plus   = rand(dist_η, Nt)
    end

    # Set n_ant_shocks shocks to 0 in pre-ZLB time periods
    if n_ant_shocks > 0
        ant1_ind = m.exogenous_shocks[:rm_shl1]
        antn_ind = m.exogenous_shocks[symbol("rm_shl$(n_ant_shocks)")]
        shock_inds = ant1_ind:antn_ind
        period_inds = vcat(inds_presample_periods(m), inds_prezlb_periods(m))
        η_plus[shock_inds, period_inds] = 0
    end

    # Produce "fake" states and observables (a+ and y+) by
    # iterating the state-space system forward
    iterate(α_plus_t1, η_plus_t) = C + T*α_plus_t1 + R*η_plus_t

    α_plus       = zeros(S, Nz, Nt)
    α_plus[:, 1] = iterate(α_plus_0, η_plus[:, 1])
    for t = 2:Nt
        α_plus[:, t] = iterate(α_plus[:, t-1], η_plus[:, t])
    end
    data_plus = D .+ Z*α_plus

    # Replace fake data with NaNs wherever actual data has NaNs
    data_plus[isnan(data)] = NaN

    # Compute y* = y - y+
    data_star = data - data_plus

    ## Run the kalman filter
    A0, P0, pred, vpred, T, R, C = if n_ant_shocks > 0

        # Note that we pass in `zeros(size(D))` instead of `D` because the
        # measurement equation for `data_star` has no constant term
        k, _, _, R3 = kalman_filter_2part(m, data_star, T, R, C, A0, P0;
            ZZ = Z, DD = zeros(size(D)), QQ = Q, MM = M, EE = E, VVall = V_all,
            allout = true, include_presample = true)

        k[:z0], k[:vz0], k[:pred], k[:vpred], R3[:TTT], R3[:RRR], R3[:CCC]
    else
        k = kalman_filter(m, data_star, T, C, Z, zeros(size(D)), V_all, A0, P0; lead = 0, allout = true)

        A0, P0, k[:pred], k[:vpred], T, R, C
    end

    ##### Step 2: Kalman smooth over everything
    α_hat_star, η_hat_star = kalman_smoother(m, data_star, T, R, C, Q, Z,
        zeros(size(D)), A0, P0, pred, vpred; include_presample = true)

    # Compute draw (states and shocks)
    α_hat = α_plus + α_hat_star
    η_hat = η_plus + η_hat_star

    if !include_presample
        α_hat = α_hat[:, index_mainsample_start(m):end]
        η_hat = η_hat[:, index_mainsample_start(m):end]
    end

    return α_hat, η_hat
end

"""
```
Hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, kal::Kalman{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

Hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, kal::Kalman{S};
    include_presample::Bool = false)

Hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, z0::Vector{S},
    pred::Matrix{S}, vpred::Array{S, 3},
    filt::Matrix{S}, vfilt::Array{S, 3};
    include_presample::Bool = false)
```
This is a Kalman Smoothing program based on the treatment in James Hamilton's
\"Time Series Analysis\". Unlike the disturbance smoother, this one does
rely on inverting singular matrices using the Moore-Penrose pseudoinverse.

Smoothed shocks are extracted by mapping the forecast errors implied by the
smoothed states back into shocks. If R is rank-deficient, then this mapping
will not exist.

### Inputs:

- `m`: model object
- `data`: the (`Ny` x `Nt`) matrix of observable data
- `T`: the (`Nz` x `Nz`) transition matrix
- `R`: the (`Nz` x `Ne`) matrix translating shocks to states
- `z0`: the (`Nz` x 1) initial (time 0) states vector
- `pred`: the (`Nz` x `Nt`) matrix of one-step-ahead predicted states (from the
  Kalman Filter)
- `vpred`: the (`Nz` x `Nz` x `Nt`) matrix of one-step-ahead predicted
  covariance matrices
- `filt`: the (`Nz` x `Nt`) matrix of filtered states
- `vfilt`: the (`Nz` x `Nz` x `Nt`) matrix of filtered covariance matrices
- `cond_type`: optional keyword argument specifying the conditional data type:
  one of `:none`, `:semi`, or `:full`. This is only necessary when a DataFrame
  (as opposed to a data matrix) is passed in, so that `df_to_matrix` knows how
  many periods of data to keep
- `include_presample`: indicates whether or not to return presample periods in
  the returned smoothed states and shocks. Defaults to `false`

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

- `α_hat`: the (`Nz` x `Nt`) matrix of smoothed states
- `η_hat`: the (`Ne` x `Nt`) matrix of smoothed shocks

If `n_presample_periods(m)` is nonzero, the `α_hat` and `η_hat` matrices will be
shorter by that number of columns (taken from the beginning).

### Notes

The state space model is defined as follows:
```
y(t) = Z*α(t) + D             (state or transition equation)
α(t+1) = T*α(t) + R*η(t+1)    (measurement or observation equation)
```
"""
function Hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, kal::Kalman{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # convert dataframe to matrix
    data = df_to_matrix(m, df, cond_type = cond_type)

    Hamilton_smoother(m, data, system, kal;
                      include_presample = include_presample)
end

function Hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, kal::Kalman{S};
    include_presample::Bool = false)

    # extract system matrices
    T, R = system[:TTT], system[:RRR]
    # extract filtered and predicted states and variances
    filt, vfilt = kal[:filt], kal[:vfilt]
    pred, vpred = kal[:pred], kal[:vpred]
    # extract initial state
    z0 = kal[:z0]

    # call actual Kalman smoother
    Hamilton_smoother(m, data, T, R, z0, pred, vpred, filt, vfilt;
        include_presample = include_presample)
end

function Hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, z0::Vector{S},
    pred::Matrix{S}, vpred::Array{S, 3},
    filt::Matrix{S}, vfilt::Array{S, 3};
    include_presample::Bool = false)

    Ne = size(R, 2)
    Ny = size(data, 1)
    Nt = size(data, 2)
    Nz = size(T, 1)

    # Check data is well-formed wrt model settings
    @assert Ny == n_observables(m)
    @assert Nt >= n_presample_periods(m) + n_prezlb_periods(m) + n_zlb_periods(m)

    # Anticipated monetary policy shocks
    n_ant_shocks = n_anticipated_shocks(m)
    t_zlb_start  = index_zlb_start(m)

    α_hat = copy(filt)

    # smooth the states recursively, starting at the Nt-1 period and going backwards
    for t in (Nt-1):-1:1
        J_t = vfilt[:,:,t] * T' * pinv(vpred[:,:,t+1])
        α_hat[:,t] = filt[:,t] + J_t * (α_hat[:,t+1] - pred[:,t+1])
    end

    # R must have rank >= n_shocks_exogenous
    if rank(R) < n_shocks_exogenous(m)
        warn("R is not sufficient rank to map forecast errors uniquely onto shocks")
    end

    η_hat = zeros((Ne,Nt))
    R_inv = pinv(R)

    # We can map the forecast errors implied by the smoothed states back to shocks
    for t in 1:Nt
        s_t = α_hat[:,t]
        s_t1 = (t == 1) ? z0 : α_hat[:,t-1]
        η_hat[:,t] = R_inv * (s_t - T * s_t1)
    end

    # trim the presample if needed
    if !include_presample
        α_hat = α_hat[:, index_mainsample_start(m):end]
        η_hat = η_hat[:, index_mainsample_start(m):end]
    end

    return α_hat, η_hat
end

# This is a Kalman Smoothing program based on the treatment in James Hamilton's
# \"Time Series Analysis\". Unlike the disturbance smoother, this one does
# rely on inverting singualr matrices using the Moore-Penrose pseudoinverse.
