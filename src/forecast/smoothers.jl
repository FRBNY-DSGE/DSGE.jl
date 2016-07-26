"""
```
function kalman_smoother{S<:AbstractFloat}(A0, P0, y, pred::Matrix{S},
    vpred::Array{S, 3}, T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Matrix{S}, n_anticipated_shocks::Int, n_anticipated_lags::Int,
    peachcount::Int, n_conditional_periods::Int, Ny0::Int = 0)
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

- `A0`: the (`Nz` x 1) initial (time 0) states vector
- `P0`: the (`Nz` x `Nz`) initial (time 0) state covariance matrix
- `y`: the (`Ny` x `Nt`) matrix of observable data
- `pred`: the (`Nz` x `Nt`) matrix of one-step-ahead predicted states (from the
  Kalman Filter)
- `vpred`: the (`Nz` x `Nz` x `Nt`) matrix of one-step-ahead predicted
  covariance matrices
- `T`: the (`Nz` x `Nz`) transition matrix
- `R`: the (`Nz` x `Ne`) matrix translating shocks to states
- `Q`: the (`Ne` x `Ne`) covariance matrix for the shocks
- `Z`: the (`Ny` x `Nz`) measurement matrix
- `D`: the (`Ny` x 1) constant vector in the measurement equation

- `n_anticipated_shocks`: an optional scalar for the zero bound specification
  indicating the number of periods ahead the interest rate is fixed
- `n_anticipated_lags`: an optional scalar for the zero bound specification
  indicating the number of periods for which interest rate expectations have
  been fixed
- `peachcount`
- `n_conditional_periods`: number of periods for which we have conditional data
  from Dick Peach
- `Ny0`: an optional scalar indicating the number of periods of presample
  (i.e. the number of periods for which smoothed states are not required)

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

`kalman_smoother` returns a `KalmanSmooth` object, which contains the following fields:
- `states`: the (`Nz` x `Nt`) matrix of smoothed states
- `shocks`: the (`Ne` x `Nt`) matrix of smoothed shocks

If `Ny0` is nonzero, the `α_hat` and `η_hat` matrices will be shorter by
that number of columns (taken from the beginning).

### Notes

The state space model is defined as follows:
```
y(t) = Z*α(t) + D             (state or transition equation) 
α(t+1) = T*α(t) + R*η(t+1)    (measurement or observation equation)
```
"""
function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Matrix{S}, A0::Array{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3})

    Ne = size(R, 2)
    Ny = size(data, 1)
    Nt = size(data, 2)
    Nz = size(T, 1)

    # Conditional data metadata
    n_unconditional_periods = subtract_quarters(date_forecast_start(m), date_presample_start(m))
    n_conditional_periods = Nt - n_unconditional_periods
    peachcount = convert(Int64, n_conditional_periods > 0)

    # Anticipated policy shocks metadata
    n_ant_shocks = n_anticipated_shocks(m)
    n_ant_lags = n_anticipated_lags(m)

    kalsmooth = disturbance_smoother(m, data, T, R, C, Q, Z, D, pred, vpred)
    r, eta_hat = kalsmooth.states, kalsmooth.shocks

    alpha_hat = zeros(Nz, Nt)
    ah_t = A0 + P0*r[:, 1]
    alpha_hat[:, 1] = ah_t

    for t = 2:Nt

        # This section relates to the zero bound framework, in which no
        # anticipated shocks are supposed to occur before the model switch.
        # In these periods, this is accomplished by setting the relevant
        # rows and columns of the Q matrix to zero. In other periods, or in
        # specifications with zero bound off (and hence with n_anticipated_shocks = 0), the
        # normal Q matrix can be used.

        if n_ant_shocks > 0
            # The first part of the conditional below pertains to the periods in
            # which zerobound is off. To specify this period, we must account
            # for (peachcount*n_conditional_periods) since peachdata is
            # augmented to y.
            # JC 11/30/10
            if t < Nt - n_ant_lags - (peachcount*n_conditional_periods)
                Q_t = zeros(Ne, Ne)
                Q_t[1:(Ne-n_ant_shocks), 1:(Ne-n_ant_shocks)] =
                    Q[1:(Ne-n_ant_shocks), 1:(Ne-n_ant_shocks)]
                ah_t = T*ah_t + R*Q_t*R'*r[:, t]
            else
                ah_t = T*ah_t + R*Q*R'*r[:, t]
            end
        else
            ah_t = T*ah_t + R*Q*R'*r[:, t]
        end

        alpha_hat[:, t] = ah_t
    end

    Nt0 = n_presample_periods(m)
    alpha_hat = alpha_hat[:, (Nt0+1):end]
    eta_hat = eta_hat[:, (Nt0+1):end]
    
    return KalmanSmooth(alpha_hat, eta_hat)
end

"""
```
function disturbance_smoother{S<:AbstractFloat}(y::Matrix{S}, pred::Matrix{S},
    vpred::Array{S,3}, T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Matrix{S}, peachcount::Int, n_conditional_periods::Int,
    n_anticipated_shocks::Int = 0, n_anticipated_lags::Int = 0)

function disturbance_smoother{S<:AbstractFloat}(y::Matrix{S}, pred::Matrix{S},
    vpred::Array{S,3}, sys::System, peachcount::Int, n_conditional_periods::Int,
    n_anticipated_shocks::Int = 0, n_anticipated_lags::Int = 0)

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

- `y`: the (`Ny` x `Nt`) matrix of observable data
- `pred`: the (`Nz` x `Nt`) matrix of one-step-ahead predicted states (from the Kalman Filter)
- `vpred`: the (`Nz` x `Nz` x `Nt`) matrix of one-step-ahead predicted covariance matrices
- `T`: the (`Nz` x `Nz`) transition matrix
- `R`: the (`Nz` x `Ne`) matrix translating shocks to states
- `Q`: the (`Ne` x `Ne`) covariance matrix for the shocks
- `Z`: the (`Ny` x `Nz`) measurement matrix
- `D`: the (`Ny` x 1) constant vector in the measurement equation
- `peachcount`
- `n_conditional_periods`: number of periods for which we have conditional data from Dick Peach

- `n_anticipated_shocks`: an optional scalar for the zero bound specification
  indicating the of periods ahead the interest rate is fixed
- `n_anticipated_lags`: an optional scalar for the zero bound specification indicating the
  number of periods for which interest rate expectations have been fixed

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

`disturbance_smoother` returns a `KalmanSmooth` object, which contains the following fields:
- `states`: the (`Nz` x `Nt`) matrix used for state smoothing
- `shocks`: the (`Ne` x `Nt`) matrix of smoothed shocks

### Notes

The state space model is defined as follows:
```
y(t) = Z*α(t) + D             (state or transition equation) 
α(t+1) = T*α(t) + R*η(t+1)    (measurement or observation equation)
```
"""
function disturbance_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Array{S}, pred::Matrix{S}, vpred::Array{S, 3})

    Nt = size(data, 2)
    Nz = size(T, 1)

    r = zeros(Nz, Nt) # holds r_{T-1}, ..., r_0
    r_t = zeros(Nz, 1)

    Ne = size(R, 2)
    eta_hat = zeros(Ne, Nt)

    # Conditional data metadata
    n_unconditional_periods = subtract_quarters(date_forecast_start(m), date_presample_start(m))
    n_conditional_periods = Nt - n_unconditional_periods
    peachcount = convert(Int64, n_conditional_periods > 0)

    # Anticipated policy shocks metadata
    n_ant_shocks = n_anticipated_shocks(m)
    n_ant_lags = n_anticipated_lags(m)

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
        # n_ant_shocks = 0), the normal Q matrix can be used.
        if n_ant_shocks > 0
            
            # The first part of the conditional below pertains to the periods in
            # which zerobound is off.  To specify this period, we must account
            # for (peachcount*n_conditional_periods) since peachdata is
            # augmented to data.  JC 11/30/10
            if t < Nt - n_ant_lags - peachcount*n_conditional_periods
                Q_t = zeros(Ne, Ne)
                Q_t[1:Ne-n_ant_shocks, 1:Ne-n_ant_shocks] =
                    Q[1:Ne-n_ant_shocks, 1:Ne-n_ant_shocks]
                eta_hat[:, t] = Q_t * R' * r_t
            else
                eta_hat[:, t] = Q * R' * r_t
            end
        else
            eta_hat[:, t] = Q * R' * r_t
        end
    end

    return KalmanSmooth(r, eta_hat)
end

function disturbance_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, sys::System, pred::Matrix{S}, vpred::Array{S, 3})

    T = sys[:TTT]
    R = sys[:RRR]
    C = sys[:CCC]
    Q  = sys[:QQ]
    Z  = sys[:ZZ]
    D  = sys[:DD]

    disturbance_smoother(m, data, T, R, C, Q, Z, D, pred, vpred)
end

"""
```
function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, P0::Matrix{S}; use_expected_rate_data = true,
    Ny0::Int = 0)

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, P0::Matrix{S}; mainsample_start = NaN, zlb_start
    = NaN, n_conditional_periods = 0, Ny0::Int = 0)
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
- `df`: the (`Nt` x `Ny`) DataFrame of observable data
- `T`: the (`Nz` x `Nz`) transition matrix
- `R`: the (`Nz` x `Ne`) matrix translating shocks to states
- `Q`: the (`Ne` x `Ne`) covariance matrix for the shocks
- `Z`: the (`Ny` x `Nz`) measurement matrix
- `D`: the (`Ny` x 1) constant vector in the measurement equation
- `P0`: the (`Nz` x `Nz`) initial (time 0) state covariance matrix. If
  `use_expected_rate_data = true`, then `P0` must include rows and columns for
  the anticipated shocks.
- `use_expected_rate_data`
- `Ny0`: an optional scalar indicating the number of periods of presample
  (i.e. the number of periods for which smoothed states are not required)

Where:

- `Nz`: number of states
- `Ny`: number of observables
- `Ne`: number of shocks
- `Nt`: number of periods for which we have data

### Outputs:

`durbin_koopman_smoother` returns a `KalmanSmooth` object, which contains the following fields:
- `states`: the (`Nz` x `Nt`) matrix of smoothed states.
- `shocks`: the (`Ne` x `Nt`) matrix of smoothed shocks.

### Notes

The state space model is defined as follows:
```
y(t) = Z*α(t) + D             (state or transition equation) 
α(t+1) = T*α(t) + R*η(t+1)    (measurement or observation equation)
```
"""
function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Matrix{S}, A0::Vector{S}, P0::Matrix{S})

    # convert DataFrame to Matrix
    data = df_to_matrix(df)
    
    # call actual simulation smoother
    durbin_koopman_smoother(m, data, T, R, C, Q, Z, D, A0, P0)
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Matrix{S}, A0::Array{S}, P0::Matrix{S})

    # Use consistent notation
    # We require the full data to be obs x periods, whereas the data frame (converted to
    # matrix) is periods x obs
    datat = data'
    Nt0 = n_presample_periods(m)
    YY0 = datat[:, 1:Nt0]
    YY  = datat[:, Nt0+1:end]

    # Get matrix dimensions
    Ny = size(YY,1) # number of observables
    Nt = size(YY,2) # number of (main-sample + ZLB) periods of data
    Nz = size(T,1)  # number of states
    Ne = size(R,2)  # number of shocks
    @assert Ny == n_observables(m)
    
    n_ant_shocks = n_anticipated_shocks(m)   # # of anticipated monetary policy shocks
    n_ant_lags   = n_anticipated_lags(m)     # # of periods for which rate expectations have been fixed
    t_zlb_start  = zlb_start_index(m)

    Nt_main      = (t_zlb_start - 1) - Nt0   # # of main-sample periods
    Nt_zlb       = Nt - Nt_main              # # of ZLB periods
   
    # Initialize matrices
    α_all_plus  = fill(NaN, Nz, Nt0+Nt)
    YY_all_plus = fill(NaN, Ny, Nt0+Nt)
    
    # Draw initial state α_0+ and sequence of shocks η+
    U, S, V = svd(P0)

    # If testing, set initial state and all shocks to zero
    if m.testing
        ap_t       = U * diagm(sqrt(S)) * zeros(Nz, 1)
        η_all_plus = sqrt(Q) * zeros(Ne, Nt0+Nt)
    else
        ap_t       = U * diagm(sqrt(S)) * randn(Nz, 1)
        η_all_plus = sqrt(Q) * randn(Ne, Nt0+Nt)
    end
    
    # Set n_ant_shocks shocks to 0 in pre-ZLB time periods
    if n_ant_shocks > 0
        # get the indices of the anticipated shocks in the m.exogenous_shocks field
        ant1_ind = m.exogenous_shocks[:rm_shl1]
        antn_ind = m.exogenous_shocks[symbol("rm_shl$(n_ant_shocks)")]

        # set shocks to 0
        η_all_plus[ant1_ind:antn_ind, 1:t_zlb_start-1] = 0
    end
    
    # Produce "fake" states and observables (a+ and y+) by
    # iterating the state-space system forward
    for t = 1:Nt0+Nt
        ap_t             = T * ap_t + R * η_all_plus[:,t]
        α_all_plus[:,t]  = ap_t
        YY_all_plus[:,t] = Z*ap_t + D
    end

    # Replace fake data with NaNs wherever actual data has NaNs
    YY_all = [YY0 YY];
    YY_all_plus[isnan(YY_all)] = NaN
    
    # Compute y* = y - y+ - D
    YY_star = YY_all - YY_all_plus

    ## Run the kalman filter
    A0, P0, pred, vpred, T, R, C, Q, Z, D = if n_ant_shocks > 0

        # PZL 2015-07-25
        # R2, R3, R1 = kalman_filter_2part(m, YY_star', allout = true, augment_states = true)
        R2, R3, R1 = kalman_filter_2part(m, YY_star', T, R, C, zeros(72, 1), P0, allout = true, augment_states = true)
        
        # unpack the results to pass to kalman_smoother

        T = R3[:TTT]
        R = R3[:RRR]
        C = R3[:CCC]

        # get expanded versions of A0 and P0
        before_shocks    = 1:(n_states(m)-n_ant_shocks)
        after_shocks_old = (n_states(m)-n_ant_shocks+1):(n_states_augmented(m)-n_ant_shocks)
        after_shocks_new = (n_states(m)+1):n_states_augmented(m)

        A0_small        = R1[:A0]
        A0              = zeros(Nz)
        A0[before_shocks]    = A0_small[before_shocks] 
        A0[after_shocks_new] = A0_small[after_shocks_old]

        P0_small        = R1[:P0]
        P0              = zeros(Nz, Nz)
        P0[before_shocks,    before_shocks]    = P0_small[before_shocks,    before_shocks]
        P0[before_shocks,    after_shocks_new] = P0_small[before_shocks,    after_shocks_old]
        P0[after_shocks_new, before_shocks]    = P0_small[after_shocks_old, before_shocks]
        P0[after_shocks_new, after_shocks_new] = P0_small[after_shocks_old, after_shocks_old]

        # concatenate all (presample, main sample, ZLB) periods together
        pred            = hcat(R1[:pred], R2[:pred], R3[:pred])
        vpred           = cat(3, R1[:vpred], R2[:vpred], R3[:vpred])

        A0, P0, pred, vpred, T, R, C, Q, Z, D
        
    else
        
        myvar = zeros(Ny+Nz,Ny+Nz)
        myvar[1:Nz,1:Nz] = R*Q*R'
        kal = kalman_filter(YY_star', 0, C, T, D, Z, myvar, A0, P0, allout=true)

        A0, P0, kal[:pred], kal[:vpred], T, R, C, Q, Z, D
        
    end

    ##### Step 2: Kalman smooth over everything
    # PZL 2016-07-25
    # kalsmooth = kalman_smoother(m, YY_star, T, R, C, Q, Z, D, A0, P0, pred, vpred)
    kalsmooth = kalman_smoother(m, YY_star, T, R, C, Q, Z, zeros(size(D)), A0, P0, pred, vpred)
    α_hat_star, η_hat_star = kalsmooth.states, kalsmooth.shocks
    
    # Compute draw (states and shocks)
    # Since `kalman_smoother` (like `durbin_koopman_smoother`) returns smoothed
    # states sans presample, we take only the main-sample and ZLB periods from
    # `α_all_plus` and `η_all_plus`
    α_til = α_all_plus[:, (Nt0+1):end] + α_hat_star
    η_til = η_all_plus[:, (Nt0+1):end] + η_hat_star

    return KalmanSmooth(α_til, η_til)
end
    
"""
```
KalmanSmooth{T<:AbstractFloat}
```

### Fields

- `states::Matrix{T}`: (`nstates` x `nperiods`) matrix of smoothed state values 
- `shocks::Matrix{T}`: (`nshocks` x `nperiods`) matrix of smoothed shock values
"""
immutable KalmanSmooth{T<:AbstractFloat}
    states::Matrix{T}
    shocks::Matrix{T}
end
