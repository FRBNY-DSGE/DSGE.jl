"""
```
kalman_smoother()
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
in the eta_hat matrix.

Nz will stand for the number of states, Ny for the number of observables,
Ne for the number of shocks, and Nt for the number of periods of data.

The state space is assumed to take the form:
y(t) = Z*alpha(t) + b
alpha(t+1) = T*alpha(t) + R*eta(t+1)

Inputs
------

- `A0`, the (Nz x 1) initial (time 0) states vector.
- `P0`, the (Nz x Nz) initial (time 0) state covariance matrix.
- `y`, the (Ny x Nt) matrix of observable data.
- `pred`, the (Nz x Nt) matrix of one-step-ahead predicted states (from the Kalman Filter).
- `vpred`, the (Nz x Nz x Nt) matrix of one-step-ahead predicted covariance matrices.
- `T`, the (Nz x Nz) transition matrix.
- `R`, the (Nz x Ne) matrix translating shocks to states.
- `Q`, the (Ne x Ne) covariance matrix for the shocks.
- `Z`, the (Ny x Nz) measurement matrix.
- `b`, the (Ny x 1) constant vector in the measurement equation.

- `n_anticipated_shocks`, an optional scalar for the zero bound specification indicating the
      number of periods ahead the interest rate is fixed.
- `antlags`, an optional scalar for the zero bound specification indicating
      the number of periods for which interest rate expectations have
      been fixed
- `Ny0`, an optional scalar indicating the number of periods of presample
      (i.e. the number of periods for which smoothed states are not required).

OUTPUTS:

- `α_hat`, the (`Nz` x `Nt`) matrix of smoothed states.
- `η_hat`, the optional (`Ne` x `Nt`) matrix of smoothed shocks.

If `Ny0` is nonzero, the `α_hat` and `η_hat` matrices will be shorter by
that number of columns (taken from the beginning).
"""
function kalman_smoother{S<:AbstractFloat}(A0, P0, y, pred::Matrix{S}, vpred::Array{S, 3}, T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S}, b::Matrix{S}, n_anticipated_shocks::Int, antlags::Int, peachcount::Int, psize::Int, Ny0::Int = 0)

    Ne = size(R, 2)
    Ny = size(y, 1)
    Nt = size(y, 2)
    Nz = size(T, 1)

    alpha_hat = zeros(Nz, Nt)

    r, eta_hat = disturbance_smoother(y, pred, vpred, T, R, Q, Z, b, peachcount, psize, n_anticipated_shocks, antlags)
    
    ah_t = A0 + P0*r[:, 1]
    alpha_hat[:, 1] = ah_t

    for t = 2:Nt

        # This section relates to the zero bound framework, in which no
        # anticipated shocks are supposed to occur before the model switch.
        # In these periods, this is accomplished by setting the relevant
        # rows and columns of the Q matrix to zero. In other periods, or in
        # specifications with zero bound off (and hence with n_anticipated_shocks = 0), the
        # normal Q matrix can be used.

        if n_anticipated_shocks != 0
            # The first part of the conditional below pertains to the periods in which zerobound is off.
            # To specify this period, we must account for (peachcount*psize) since peachdata is augmented to y.
            # JC 11/30/10
            if n_anticipated_shocks > 0 && t < Nt - antlags - (peachcount*psize)
                Q_t = zeros(Ne, Ne)
                Q_t[1:(Ne-n_anticipated_shocks), 1:(Ne-n_anticipated_shocks)] = Q[1:(Ne-n_anticipated_shocks), 1:(Ne-n_anticipated_shocks)]
                ah_t = T*ah_t + R*Q_t*R'*r[:, t]
            else
                ah_t = T*ah_t + R*Q*R'*r[:, t]
            end
        else
            ah_t = T*ah_t + R*Q*R'*r[:, t]
        end

        alpha_hat[:, t] = ah_t
    end

    alpha_hat = alpha_hat[:, (Ny0+1):end]
    eta_hat = eta_hat[:, (Ny0+1):end]
    
    return alpha_hat, eta_hat
end





# DISTSMTH_K93.M

# This is a Kalman Smoothing program based on S.J. Koopman's "Disturbance
# Smoother for State Space Models" (Biometrika, 1993), as specified in
# Durbin and Koopman's "A Simple and Efficient Simulation Smoother for
# State Space Time Series Analysis" (Biometrika, 2002). The algorithm has been
# simplified for the case in which there is no measurement error, and the
# model matrices do not vary with time.

# This disturbance smoother is intended for use with the state smoother
# kalsmth_93.m from the same papers (Koopman 1993, Durbin and Koopman
# 2002). It produces a matrix of vectors, r, that is used for state
# smoothing, and an optional matrix, eta_hat, containing the smoothed
# shocks. It has been adjusted to account for the possibility of missing
# values in the data, and to accommodate the zero bound model, which
# requires that no anticipated shocks occur before the zero bound window,
# which is achieved by setting the entries in the Q matrix corresponding to
# the anticipated shocks to zero in those periods.

# Nz will stand for the number of states, Ny for the number of observables,
# Ne for the number of shocks, and Nt for the number of periods of data.

# The state space is assumed to take the form:
# y(t) = Z*alpha(t) + b
# alpha(t+1) = T*alpha(t) + R*eta(t+1)

# INPUTS:

# y, the (Ny x Nt) matrix of observable data.
# pred, the (Nz x Nt) matrix of one-step-ahead predicted states (from the Kalman Filter).
# vpred, the (Nz x Nz x Nt) matrix of one-step-ahead predicted covariance matrices.
# T, the (Nz x Nz) transition matrix.
# R, the (Nz x Ne) matrix translating shocks to states.
# Q, the (Ne x Ne) covariance matrix for the shocks.
# Z, the (Ny x Nz) measurement matrix.
# b, the (Ny x 1) constant vector in the measurement equation.

# n_anticipated_shocks, an optional scalar for the zero bound specification indicating the
#       number of periods ahead the interest rate is fixed.
# antlags, an optional scalar for the zero bound specification indicating
#       the number of periods for which interest rate expectations have
#       been fixed
# Ny0, an optional scalar indicating the number of periods of presample
#       (i.e. the number of periods for which smoothed states are not
#       required).

# OUTPUTS:

# r, the (Nz x Nt) matrix used for state smoothing.
# eta_hat, the optional (Ne x Nt) matrix of smoothed shocks.

# Dan Greenwald, 7/7/2010.
function disturbance_smoother{S<:AbstractFloat}(y::Matrix{S}, pred::Matrix{S}, vpred::Array{S,3}, T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S}, b::Matrix{S}, peachcount::Int, psize::Int, n_anticipated_shocks::Int = 0, antlags::Int = 0)
    Nt = size(y, 2)
    Nz = size(T, 1)

    r = zeros(Nz, Nt) # holds r_T-1, ...r_0
    r_t = zeros(Nz, 1)

    Ne = size(R, 2)
    eta_hat = zeros(Ne, Nt)

    for t = Nt:-1:1

        y_t = y[:, t]

        # This section deals with the possibility of missing values in the y_t
        # vector (especially relevant for smoothing over peachdata).
        nonmissing = !isnan(y_t)
        y_t = y_t[nonmissing]
        Z_t = Z[nonmissing, :]
        b_t = b[nonmissing]

        a = pred[:, t]
        P = vpred[:, :, t]

        F = Z_t*P*Z_t'
        v = y_t - Z_t*a - b_t
        K = T*P*Z_t'/F
        L = T - K*Z_t

        r_t = Z_t'/F*v + L'*r_t
        r[:, t] = r_t

        # This section relates to the zero bound framework, in which no
        # anticipated shocks are supposed to occur before the model switch.
        # In these periods, this is accomplished by setting the relevant
        # rows and columns of the Q matrix to zero. In other periods, or in
        # specifications with zero bound off (and hence with n_anticipated_shocks = 0), the
        # normal Q matrix can be used.
        if n_anticipated_shocks != 0
            
            # The first part of the conditional below pertains to the periods in which zerobound is off.
            # To specify this period, we must account for (peachcount*psize) since peachdata is augmented to y.
            # JC 11/30/10
            if n_anticipated_shocks > 0 && t < Nt - antlags - peachcount*psize
                Q_t = zeros(Ne, Ne)
                Q_t[1:Ne-n_anticipated_shocks, 1:Ne-n_anticipated_shocks] = Q[1:Ne-n_anticipated_shocks, 1:Ne-n_anticipated_shocks]
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

function disturbance_smoother{S<:AbstractFloat}(y::Matrix{S}, pred::Matrix{S}, vpred::Array{S,3},
                                                sys::System, peachcount::Int, psize::Int,
                                                n_anticipated_shocks::Int = 0, antlags::Int = 0)

    TTT = sys[:TTT]
    RRR = sys[:RRR]
    QQ  = sys[:QQ]
    ZZ  = sys[:ZZ]
    DD  = sys[:DD]

    disturbance_smoother(y, pred, vpred, TTT, RRR, QQ, ZZ, DD, peachcount, psize, n_anticipated_shocks, antlags)
end


    
"""
```
KalmanSmooth{T<:AbstractFloat}

Fields
------
- `states::Matrix{T}`: nstates x nperiods matrix of smoothed state values 
- `shocks::Matrix{T}`: nshocks x nperiods matrix of smoothed shock values
```
"""
immutable KalmanSmooth{T<:AbstractFloat}
    states::Matrix{T}
    shocks::Matrix{T}
end


