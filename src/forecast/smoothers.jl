"""
```
kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, z0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    cond_type::Symbol = :none, include_presample::Bool = false)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, z0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    include_presample::Bool = false)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, z0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    cond_type::Symbol = :none, include_presample::Bool = false)

kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S}
    T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S}, Z::Matrix{S},
    D::Vector{S}, z0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
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
- `z0`: the (`Nz` x 1) initial (time 0) states vector
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
function kalman_smoother{S<:AbstractFloat}(data::Matrix{S},
    TTT::Matrix{S}, RRR::Matrix{S}, CCC::Vector{S},
    QQ::Matrix{S}, ZZ::Matrix{S}, DD::Vector{S},
    z0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_presample_periods::Int = 0)

    T = size(data, 2)
    regime_indices = Range{Int64}[1:T]

    kalman_smoother(regime_indices, data, Matrix{S}[TTT], Matrix{S}[RRR], Vector{S}[CCC],
        Matrix{S}[QQ], Matrix{S}[ZZ], Vector{S}[DD], z0, P0, pred, vpred;
        n_presample_periods = n_presample_periods)
end

function kalman_smoother{S<:AbstractFloat}(regime_indices::Vector{Range{Int64}},
    data::Matrix{S}, TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}},
    QQs::Vector{Matrix{S}}, ZZs::Vector{Matrix{S}}, DDs::Vector{Vector{S}},
    z0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    n_presample_periods::Int = 0)

    n_regimes = length(regime_indices)

    # Dimensions
    T  = size(data,    2) # number of periods of data
    Nz = size(TTTs[1], 1) # number of states

    # Call disturbance smoother
    smoothed_disturbances, smoothed_shocks = disturbance_smoother(regime_indices, data,
                                                 TTTs, RRRs, QQs, ZZs, DDs, pred, vpred)

    # Initialize outputs
    smoothed_states = zeros(S, Nz, T)

    r     = smoothed_disturbances[:, 1]
    α_hat = z0 + P0*r
    smoothed_states[:, 1] = α_hat

    for i = 1:n_regimes
        # Get state-space system matrices for this regime
        regime_periods = regime_indices[i]

        TTT, RRR     = TTTs[i], RRRs[i]
        QQ,  ZZ,  DD = QQs[i],  ZZs[i],  DDs[i]

        for t in regime_periods
            # t = 1 has already been initialized
            t == 1 ? continue : nothing

            r     = smoothed_disturbances[:, t]
            α_hat = TTT*α_hat + RRR*QQ*RRR'*r

            smoothed_states[:, t] = α_hat
        end
    end

    if n_presample_periods > 0
        mainsample_periods = n_presample_periods+1:T

        smoothed_states = smoothed_states[:, mainsample_periods]
        smoothed_shocks = smoothed_shocks[:, mainsample_periods]
    end

    return smoothed_states, smoothed_shocks
end

function kalman_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, z0::Vector{S}, P0::Matrix{S}, pred::Matrix{S}, vpred::Array{S, 3};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # Partition sample into pre- and post-ZLB regimes
    # Note that the post-ZLB regime may be empty if we do not impose the ZLB
    regime_inds = zlb_regime_indices(m, data)

    # Get system matrices for each regime
    TTTs, RRRs, CCCs, QQs, ZZs, DDs, _ = zlb_regime_matrices(m, system)

    # Specify number of presample periods if we don't want to include them in
    # the final results
    T0 = include_presample ? 0 : n_presample_periods(m)

    # Call Kalman smoother
    kalman_smoother(regime_inds, data, TTTs, RRRs, CCCs,
        QQs, ZZs, DDs, z0, P0, pred, vpred;
        n_presample_periods = T0)
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
function disturbance_smoother{S<:AbstractFloat}(data::Matrix{S},
    TTT::Matrix{S}, RRR::Matrix{S},
    QQ::Matrix{S}, ZZ::Matrix{S}, DD::Vector{S},
    pred::Matrix{S}, vpred::Array{S, 3},
    n_presample_periods::Int = 0)

    T = size(data, 2)
    regime_indices = Range{Int64}[1:T]

    disturbance_smoother(regime_indices, data, Matrix{S}[TTT], Matrix{S}[RRR],
        Matrix{S}[QQ], Matrix{S}[ZZ], Vector{S}[DD], pred, vpred;
        n_presample_periods = n_presample_periods)
end

function disturbance_smoother{S<:AbstractFloat}(regime_indices::Vector{Range{Int64}},
    data::Matrix{S}, TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}},
    QQs::Vector{Matrix{S}}, ZZs::Vector{Matrix{S}}, DDs::Vector{Vector{S}},
    pred::Matrix{S}, vpred::Array{S, 3}; n_presample_periods::Int = 0)

    n_regimes = length(regime_indices)

    # Dimensions
    T  = size(data,    2) # number of periods of data
    Nz = size(TTTs[1], 1) # number of states
    Ne = size(RRRs[1], 2) # number of shocks

    # Initialize outputs
    smoothed_disturbances = zeros(S, Nz, T)
    smoothed_shocks       = zeros(S, Ne, T)

    r = zeros(S, Nz)

    for i = n_regimes:-1:1
        # Get state-space system matrices for this regime
        regime_periods = regime_indices[i]

        TTT, RRR     = TTTs[i], RRRs[i]
        QQ,  ZZ,  DD = QQs[i],  ZZs[i],  DDs[i]

        for t in reverse(regime_periods)
            # If an element of the vector y_t is missing (NaN) for the observation t, the
            # corresponding row is ditched from the measurement equation
            nonmissing = !isnan(data[:, t])
            y_t  = data[nonmissing, t]
            ZZ_t = ZZ[nonmissing, :]
            DD_t = DD[nonmissing]

            a = pred[:, t]
            P = vpred[:, :, t]

            F = ZZ_t*P*ZZ_t'
            v = y_t - ZZ_t*a - DD_t
            K = TTT*P*ZZ_t'/F
            L = TTT - K*ZZ_t

            r = ZZ_t'/F*v + L'*r
            smoothed_disturbances[:, t] = r

            smoothed_shocks[:, t] = QQ*RRR'*r

        end # of loop backward through this regime's periods

    end # of loop backward through regimes

    if n_presample_periods > 0
        mainsample_periods = n_presample_periods+1:T

        smoothed_disturbances = smoothed_disturbances[:, mainsample_periods]
        smoothed_shocks       = smoothed_shocks[:,       mainsample_periods]
    end

    return smoothed_disturbances, smoothed_shocks
end

"""
```
durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, system::System, z0::Vector{S}, P0::Matrix{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data:Matrix{S}, system::System, z0::Vector{S}, P0::Matrix{S};
    include_presample::Bool = false)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    df::DataFrame, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, z0::Array{S}, P0::Matrix{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Array{S}, Q::Matrix{S},
    Z::Matrix{S}, D::Matrix{S}, z0::Array{S}, P0::Matrix{S};
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
- `z0`: the (`Nz` x 1) initial (time 0) states vector
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
function durbin_koopman_smoother{S<:AbstractFloat}(data::Matrix{S},
    TTT::Matrix{S}, RRR::Matrix{S}, CCC::Vector{S},
    QQ::Matrix{S}, ZZ::Matrix{S}, DD::Vector{S},
    MM::Matrix{S}, EE::Matrix{S}, z0::Vector{S}, P0::Matrix{S};
    n_presample_periods::Int = 0, draw_states::Bool = true)

    T = size(data, 2)
    regime_indices = Range{Int64}[1:T]

    durbin_koopman_smoother(regime_indices, data, Matrix{S}[TTT], Matrix{S}[RRR], Vector{S}[CCC],
        Matrix{S}[QQ], Matrix{S}[ZZ], Vector{S}[DD], Vector{S}[MM], Vector{S}[EE], z0, P0;
        n_presample_periods = n_presample_periods, draw_states = draw_states)
end

function durbin_koopman_smoother{S<:AbstractFloat}(regime_indices::Vector{Range{Int64}},
    data::Matrix{S}, TTTs::Vector{Matrix{S}}, RRRs::Vector{Matrix{S}}, CCCs::Vector{Vector{S}},
    QQs::Vector{Matrix{S}}, ZZs::Vector{Matrix{S}}, DDs::Vector{Vector{S}},
    MMs::Vector{Matrix{S}}, EEs::Vector{Matrix{S}}, z0::Vector{S}, P0::Matrix{S};
    n_presample_periods::Int = 0, draw_states::Bool = true)

    n_regimes = length(regime_indices)

    # Dimensions
    T  = size(data,    2) # number of periods of data
    Nz = size(TTTs[1], 1) # number of states
    Ne = size(RRRs[1], 2) # number of shocks
    Ny = size(ZZs[1],  1) # number of observables

    # Draw initial state α_0+ and sequence of shocks η+
    if draw_states
        U, eig, _ = svd(P0)
        α_plus_t  = U * diagm(sqrt(eig)) * randn(Nz)
        η_plus    = sqrt(QQs[1]) * randn(Ne, Nt)
    else
        α_plus_t  = zeros(S, Nz)
        η_plus    = zeros(S, Ne, T)
    end

    # Produce "fake" states and observables (α+ and y+) by
    # iterating the state-space system forward
    α_plus       = zeros(S, Nz, T)
    y_plus       = zeros(S, Ny, T)

    for i = 1:n_regimes
        # Get state-space system matrices for this regime
        regime_periods = regime_indices[i]

        TTT, RRR, CCC = TTTs[i], RRRs[i], CCCs[i]
        QQ,  ZZ,  DD  = QQs[i],  ZZs[i],  DDs[i]

        for t in regime_periods
            η_plus_t = η_plus[:, t]
            α_plus_t = TTT*α_plus_t + RRR*η_plus_t + CCC

            α_plus[:, t] = α_plus_t
            y_plus[:, t] = ZZ*α_plus_t + DD
        end
    end

    # Replace fake data with NaNs wherever actual data has NaNs
    y_plus[isnan(data)] = NaN

    # Compute y* = y - y+
    y_star = data - y_plus

    # Run the Kalman filter
    # Note that we pass in `zeros(size(D))` instead of `D` because the
    # measurement equation for `data_star` has no constant term
    _, _, _, pred, vpred, _ = kalman_filter(regime_indices, y_star, TTTs, RRRs, CCCs,
                                  QQs, ZZs, fill(zeros(Ny), n_regimes), MMs, EEs,
                                  z0, P0; allout = true)

    # Kalman smooth
    α_hat_star, η_hat_star = kalman_smoother(regime_indices, y_star, TTTs, RRRs, CCCs,
                                 QQs, ZZs, fill(zeros(Ny), n_regimes),
                                 z0, P0, pred, vpred)

    # Compute draw (states and shocks)
    smoothed_states = α_plus + α_hat_star
    smoothed_shocks = η_plus + η_hat_star

    if n_presample_periods > 0
        mainsample_periods = n_presample_periods+1:T

        smoothed_states = smoothed_states[:, mainsample_periods]
        smoothed_shocks = smoothed_shocks[:, mainsample_periods]
    end

    return smoothed_states, smoothed_shocks
end

function durbin_koopman_smoother{S<:AbstractFloat}(m::AbstractModel,
    data::Matrix{S}, system::System, z0::Vector{S}, P0::Matrix{S};
    include_presample::Bool = false)

    # Partition sample into pre- and post-ZLB regimes
    # Note that the post-ZLB regime may be empty if we do not impose the ZLB
    regime_inds = zlb_regime_indices(m, data)

    # Get system matrices for each regime
    TTTs, RRRs, CCCs, QQs, ZZs, DDs, MMs, EEs = zlb_regime_matrices(m, system)

    # Specify number of presample periods if we don't want to include them in
    # the final results
    T0 = include_presample ? 0 : n_presample_periods(m)

    # Call Durbin-Koopman smoother
    durbin_koopman_smoother(regime_inds, data, TTTs, RRRs, CCCs,
        QQs, ZZs, DDs, MMs, EEs, z0, P0;
        n_presample_periods = T0, draw_states = !m.testing)
end

"""
```
hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, kal::Kalman{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, kal::Kalman{S};
    include_presample::Bool = false)

hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, z0::Vector{S},
    pred::Matrix{S}, vpred::Array{S, 3},
    filt::Matrix{S}, vfilt::Array{S, 3};
    include_presample::Bool = false)
```
This is a Kalman Smoothing program based on the treatment in James Hamilton's
\"Time Series Analysis\". Unlike the disturbance smoother, this one does
rely on inverting potentially singular matrices using the Moore-Penrose pseudoinverse.

Smoothed shocks are extracted by mapping the forecast errors implied by the
smoothed states back into shocks. As such, this routine assumes that R has
sufficient rank to have a left inverse (i.e. that there are more states than shocks).

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
function hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, kal::Kalman{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # convert dataframe to matrix
    data = df_to_matrix(m, df, cond_type = cond_type)

    hamilton_smoother(m, data, system, kal;
                      include_presample = include_presample)
end

function hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
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
    hamilton_smoother(m, data, T, R, z0, pred, vpred, filt, vfilt;
        include_presample = include_presample)
end

function hamilton_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
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

    # smooth the states recursively, starting at the Nt-1 period and going backwards
    α_hat = copy(filt)

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

"""
```
carter_kohn_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, kal::Kalman{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

carter_kohn_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    system::System, kal::Kalman{S};
    include_presample::Bool = false)

carter_kohn_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, z0::Vector{S},
    pred::Matrix{S}, vpred::Array{S, 3},
    filt::Matrix{S}, vfilt::Array{S, 3};
    include_presample::Bool = false)
```
This program is a simulation smoother based on Carter and Kohn's
\"On Gibbs Sampling for State Space Modeks\" (Biometrika, 1994).
It recursively sampling from the conditional distribution of time t
states given the full set of observables and states from time t+1 to
time T. Unlike the Durbin Koopman simulation smoother, this one does
rely on inverting potentially singular matrices using the Moore-Penrose
pseudoinverse.

Smoothed shocks are extracted by mapping the forecast errors implied by the
smoothed states back into shocks. As such, this routine assumes that R has
sufficient rank to have a left inverse (i.e. that there are more states than shocks).

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
function carter_kohn_smoother{S<:AbstractFloat}(m::AbstractModel, df::DataFrame,
    system::System, kal::Kalman{S};
    cond_type::Symbol = :none, include_presample::Bool = false)

    # convert dataframe to matrix
    data = df_to_matrix(m, df, cond_type = cond_type)

    carter_kohn_smoother(m, data, system, kal;
                      include_presample = include_presample)
end

function carter_kohn_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
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
    carter_kohn_smoother(m, data, T, R, z0, pred, vpred, filt, vfilt;
        include_presample = include_presample)
end

function carter_kohn_smoother{S<:AbstractFloat}(m::AbstractModel, data::Matrix{S},
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

    # draw from distribution of α_t|T = α_T|T
    U, D, _ = svd(vfilt[:,:,Nt])
    dist_α = DegenerateMvNormal(filt[:,Nt], U * diagm(sqrt(D)))
    α_TT = m.testing ? filt[:,Nt] : rand(dist_α)

    # smooth the states recursively, starting at the Nt-1 period and going backwards
    α_hat = copy(filt)
    α_hat[:,Nt] = α_TT
    for t in (Nt-1):-1:1
        J_t = vfilt[:,:,t] * T' * pinv(vpred[:,:,t+1])
        μ = filt[:,t] + J_t * (α_hat[:,t+1] - pred[:,t+1])
        Σ = vfilt[:,:,t] - J_t * T * vfilt[:,:,t]

        U, D, _ = svd(Σ)
        dist_α = DegenerateMvNormal(μ, U * diagm(sqrt(D)))
        α_hat[:,t] = m.testing ? μ : rand(dist_α)
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