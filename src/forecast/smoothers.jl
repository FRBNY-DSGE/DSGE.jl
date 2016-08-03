"""
```
kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    sys::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    sys::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)
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
- `n_conditional_periods`: optional argument indicating the number of periods of
  conditional data in `data`

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

- `alpha_hat`: the (`Nz` x `Nt`) matrix of smoothed states
- `eta_hat`: the (`Ne` x `Nt`) matrix of smoothed shocks

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
    sys::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)

    # extract system matrices
    T, R, C = sys[:TTT], sys[:RRR], sys[:CCC]
    Q, Z, D = sys[:QQ], sys[:ZZ], sys[:DD]
    
    # call actual Kalman smoother
    kalman_smoother(m, df, T, R, C, Q, Z, D, A0, P0, pred, vpred;
        n_conditional_periods = n_conditional_periods)
end

function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    sys::System, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)

    # extract system matrices
    T, R, C = sys[:TTT], sys[:RRR], sys[:CCC]
    Q, Z, D = sys[:QQ], sys[:ZZ], sys[:DD]
    
    # call actual Kalman smoother
    kalman_smoother(m, data, T, R, C, Q, Z, D, A0, P0, pred, vpred;
        n_conditional_periods = n_conditional_periods)
end

function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)

    # convert DataFrame to matrix
    data = df_to_matrix(df)'
    
    # call actual Kalman smoother
    kalman_smoother(m, data, T, R, C, Q, Z, D, A0, P0, pred, vpred;
        n_conditional_periods = n_conditional_periods)
end

function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, A0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_conditional_periods::Int = 0)

    Ne = size(R, 2)
    Ny = size(data, 1)
    Nt = size(data, 2)
    Nz = size(T, 1)
    
    # Check data is well-formed wrt model settings
    @assert Ny == n_observables(m)
    @assert Nt == n_presample_periods(m) + n_prezlb_periods(m) + n_zlb_periods(m) + n_conditional_periods

    # Anticipated monetary policy shocks
    n_ant_shocks = n_anticipated_shocks(m)
    t_zlb_start  = index_zlb_start(m)

    r, eta_hat = disturbance_smoother(m, data, T, R, C, Q, Z, D, pred, vpred)

    alpha_hat = zeros(Nz, Nt)
    ah_t = A0 + P0*r[:, 1]
    alpha_hat[:, 1] = ah_t

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

        alpha_hat[:, t] = ah_t
    end

    period_inds = [inds_prezlb_periods(m); inds_zlb_periods(m)]
    alpha_hat = alpha_hat[:, period_inds]
    eta_hat   = eta_hat[:,   period_inds]
    
    return alpha_hat, eta_hat
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
- `eta_hat`: the (`Ne` x `Nt`) matrix of smoothed shocks

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
    eta_hat = zeros(Ne, Nt)

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
                eta_hat[:, t] = Q_t * R' * r_t
            else
                eta_hat[:, t] = Q * R' * r_t
            end
        else
            eta_hat[:, t] = Q * R' * r_t
        end
    end

    return r, eta_hat
end

"""
```
durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, sys::System, A0::Vector{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data:Matrix{S}, sys::System, A0::Vector{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, A0::Array{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, A0::Array{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)
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
- `n_conditional_periods`: optional argument indicating the number of periods of
  conditional data in `data`

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

- `alpha_hat`: the (`Nz` x `Nt`) matrix of smoothed states.
- `eta_hat`: the (`Ne` x `Nt`) matrix of smoothed shocks.

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
    df::DataFrame, sys::System, A0::Vector{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)

    # extract system matrices
    T, R, C = sys[:TTT], sys[:RRR], sys[:CCC]
    Q, Z, D = sys[:QQ], sys[:ZZ], sys[:DD]
    
    # call actual Durbin-Koopman smoother
    durbin_koopman_smoother(m, df, T, R, C, Q, Z, D, A0, P0;
        n_conditional_periods = n_conditional_periods)
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, sys::System, A0::Vector{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)

    # extract system matrices
    T, R, C = sys[:TTT], sys[:RRR], sys[:CCC]
    Q, Z, D = sys[:QQ], sys[:ZZ], sys[:DD]
    
    # call actual Durbin-Koopman smoother
    durbin_koopman_smoother(m, data, T, R, C, Q, Z, D, A0, P0;
        n_conditional_periods = n_conditional_periods)
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Vector{S}, A0::Vector{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)

    # convert DataFrame to Matrix
    data = df_to_matrix(df)'
    
    # call actual simulation smoother
    durbin_koopman_smoother(m, data, T, R, C, Q, Z, D, A0, P0;
        n_conditional_periods = n_conditional_periods)
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Vector{S}, A0::Array{S}, P0::Matrix{S};
    n_conditional_periods::Int = 0)

    # Get matrix dimensions
    Ny = size(data, 1)
    Nt = size(data, 2)
    Nz = size(T, 1)
    Ne = size(R, 2)
    
    # Check data is well-formed wrt model settings
    @assert Ny == n_observables(m)
    @assert Nt == n_presample_periods(m) + n_prezlb_periods(m) + n_zlb_periods(m) + n_conditional_periods

    # Anticipated monetary policy shocks
    n_ant_shocks = n_anticipated_shocks(m)
    t_zlb_start  = index_zlb_start(m)
   
    # Initialize matrices
    α_all_plus  = fill(NaN, Nz, Nt)
    YY_all_plus = fill(NaN, Ny, Nt)
    
    # Draw initial state α_0+ and sequence of shocks η+
    U, E, V = svd(P0)

    # If testing, set initial state and all shocks to zero
    if m.testing
        ap_t       = U * diagm(sqrt(E)) * zeros(Nz, 1)
        η_all_plus = sqrt(Q) * zeros(Ne, Nt)
    else
        ap_t       = U * diagm(sqrt(E)) * randn(Nz, 1)
        η_all_plus = sqrt(Q) * randn(Ne, Nt)
    end
    
    # Set n_ant_shocks shocks to 0 in pre-ZLB time periods
    if n_ant_shocks > 0
        # get the indices of the anticipated shocks in the m.exogenous_shocks
        # field
        ant1_ind = m.exogenous_shocks[:rm_shl1]
        antn_ind = m.exogenous_shocks[symbol("rm_shl$(n_ant_shocks)")]
        shock_inds = ant1_ind:antn_ind
        period_inds = [inds_presample_periods(m); inds_prezlb_periods(m)]

        # set shocks to 0
        η_all_plus[shock_inds, period_inds] = 0
    end
    
    # Produce "fake" states and observables (a+ and y+) by
    # iterating the state-space system forward
    for t = 1:Nt
        ap_t             = T * ap_t + R * η_all_plus[:,t]
        α_all_plus[:,t]  = ap_t
        YY_all_plus[:,t] = Z*ap_t + D
    end

    # Replace fake data with NaNs wherever actual data has NaNs
    YY_all_plus[isnan(data)] = NaN
    
    # Compute y* = y - y+ - D
    YY_star = data - YY_all_plus

    ## Run the kalman filter
    A0, P0, pred, vpred, T, R, C = if n_ant_shocks > 0

        # Note that we pass in `zeros(size(D))` instead of `D` because the
        # measurement equation for `YY_star` has no constant term
        k, _, _, R3 = kalman_filter_2part(m, YY_star', T, R, C, A0, P0;
            DD = zeros(size(D)), allout = true, include_presample = true)
        
        k[:z0], k[:vz0], k[:pred], k[:vpred], R3[:TTT], R3[:RRR], R3[:CCC]
    else
        VVall = zeros(Ny+Nz,Ny+Nz)
        VVall[1:Nz,1:Nz] = R*Q*R'
        
        k = kalman_filter(m, YY_star, T, C, Z, D, VVall, A0, P0; lead = 0, allout = true)

        A0, P0, k[:pred], k[:vpred], T, R, C
    end

    ##### Step 2: Kalman smooth over everything
    α_hat_star, η_hat_star = kalman_smoother(m, YY_star, T, R, C, Q, Z,
        zeros(size(D)), A0, P0, pred, vpred)
    
    # Compute draw (states and shocks)
    alpha_hat = α_all_plus[:, index_prezlb_start(m):end] + α_hat_star
    eta_hat   = η_all_plus[:, index_prezlb_start(m):end] + η_hat_star

    return alpha_hat, eta_hat
end
