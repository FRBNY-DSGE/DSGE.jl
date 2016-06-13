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
function kalman_smoother{S<:AbstractFloat}(A0, P0, y, pred::Matrix{S}, vpred::Array{S, 3},
    T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S}, b::Matrix{S},
    n_anticipated_shocks::Int, antlags::Int, peachcount::Int, psize::Int, Ny0::Int = 0)

    Ne = size(R, 2)
    Ny = size(y, 1)
    Nt = size(y, 2)
    Nz = size(T, 1)

    alpha_hat = zeros(Nz, Nt)

    r, eta_hat = disturbance_smoother(y, pred, vpred, T, R, Q, Z, b, peachcount, psize,
        n_anticipated_shocks, antlags)
    
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
function disturbance_smoother_k93{S<:AbstractFloat}(y::Matrix{S}, pred::Matrix{S},
    vpred::Array{S,3}, T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S}, b::Matrix{S},
    peachcount::Int, psize::Int, n_anticipated_shocks::Int = 0, antlags::Int = 0)
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

function disturbance_smoother{S<:AbstractFloat}(y::Matrix{S}, pred::Matrix{S},
    vpred::Array{S,3}, sys::System, peachcount::Int, psize::Int, n_anticipated_shocks::Int =
    0, antlags::Int = 0)

    TTT = sys[:TTT]
    RRR = sys[:RRR]
    QQ  = sys[:QQ]
    ZZ  = sys[:ZZ]
    DD  = sys[:DD]

    disturbance_smoother(y, pred, vpred, TTT, RRR, QQ, ZZ, DD, peachcount, psize, n_anticipated_shocks, antlags)
end


# DRAWSTATES_DK02.M

# This program is a simulation smoother based on Durbin and Koopman's
# "A Simple and Efficient Simulation Smoother for State Space Time Series
# Analysis" (Biometrika, 2002). The algorithm has been simplified for the
# case in which there is no measurement error, and the model matrices do
# not vary with time.
    
# Unlike other simulation smoothers (for example, that of Carter and Kohn,
# 1994), this method does not require separate draws for each period, draws
# of the state vectors, or even draws from a conditional distribution.
# Instead, vectors of shocks are drawn from the unconditional distribution
# of shocks, which is then corrected (via a Kalman Smoothing step), to
# yield a draw of shocks conditional on the data. This is then used to
# generate a draw of states conditional on the data. Drawing the states in
# this way is much more efficient than other methods, as it avoids the need
# for multiple draws of state vectors (requiring singular value
# decompositions), as well as inverting state covariance matrices
# (requiring the use of the computationally intensive and relatively
# erratic Moore-Penrose pseudoinverse).

# Nz will stand for the number of states, Ny for the number of observables,
# Ne for the number of shocks, and Nt for the number of periods of data.
        
# The state space is assumed to take the form:
# y(t) = Z*alpha(t) + b
# alpha(t+1) = T*alpha(t) + R*eta(t+1)

# INPUTS:

# df, the (Nt x Ny) DataFrame of observable data.
# T, the (Nz x Nz) transition matrix.
# R, the (Nz x Ne) matrix translating shocks to states.
# Q, the (Ne x Ne) covariance matrix for the shocks.
# Z, the (Ny x Nz) measurement matrix.
# b, the (Ny x 1) constant vector in the measurement equation.
# A0, the (Nz x 1) initial (time 0) states vector.
# P0, the (Nz x Nz) initial (time 0) state covariance matrix.

# nant, an optional scalar for the zero bound specification indicating the
#       number of periods ahead the interest rate is fixed.
# antlags, an optional scalar for the zero bound specification indicating
#       the number of periods for which interest rate expectations have
#       been fixed
# Ny0, an optional scalar indicating the number of periods of presample
#       (i.e. the number of periods for which smoothed states are not
#       required). If you want to set a value for Ny0 but not for nant or
#       antlags then set them both empty.                 
function drawstates_dk02!{T<:AbstractFloat}(m::AbstractModel,
                                            df::DataFrame,
                                            TTT::Matrix{T},
                                            RRR::Matrix{T},
                                            CCC::Array{T},
                                            QQ::Matrix{T},
                                            ZZ::Matrix{T},
                                            DD::Matrix{T},
                                            A0::Vector{T},
                                            P0::Matrix{T})

    ## extract settings/dates/etc

    n_ant_shocks  = n_anticipated_shocks(m)
    zlb_start_ind = zlb_start_index(m)
    mainsample_start_ind = Nt0 + 1

    
    #conditional_start_ind =
    t0 = first_forecast_quarter(m)
    t1 = df[end, :date]
    n_conditional_periods = subtract_quarters(t1, t0)   
    
    # convert DataFrame to Matrix
    data = df_to_matrix(df)
    
    # call actual simulation smoother
    drawstates_dk02!(m, data, A0, P0,
                    mainsample_start = mainsample_start_ind,
                    zlb_start = zlb_start_ind,
                    n_conditional_periods = n_conditional_periods)
end
 

function drawstates_dk02!{T<:AbstractFloat}(m::AbstractModel,
                                            data::Matrix{T},
                                            TTT::Matrix{T},
                                            RRR::Matrix{T},
                                            CCC::Array{T},
                                            QQ::Matrix{T},
                                            ZZ::Matrix{T},
                                            DD::Matrix{T},
                                            A0::Matrix{T},
                                            P0::Matrix{T};
                                            mainsample_start = NaN,
                                            zlb_start = NaN,
                                            n_conditional_periods = 0)
    
    # Use consistent notation
    # We require the full data to be obs x periods, whereas the data frame (converted to
    # matrix) is periods x obs
    datat = data'
    Nt0 = n_presample_periods(m)
    YY0 = datat[:, 1:Nt0]
    YY  = datat[:, Nt0+1:end]

    # Get matrix dimensions
    
    Ny = size(YY,1)        # number of observables
    Nt = size(YY,2)        # number of periods of data (minus presample)
    Nz = size(TTT,1)         # number of states
    Ne = size(RRR,2)         # number of shocks
    n_ant_shocks = n_anticipated_shocks(m)   # # of anticipated monetary policy shocks

    println("data: $(size(datat))")
    println("Ny: $(Ny)")
    println("Nt0: $(Nt0)")
    println("Nt: $(Nt)")
    println("Nz: $(Nz)")
    println("Ne: $(Ne)")

   
    # intialize matrices
    α_all_plus = fill(NaN, Nz, Nt0+Nt)
    YY_all_plus  = fill(NaN, Ny, Nt0+Nt)
    
    # Draw initial state, a₀+
    U,D,V = svd(P0)
    println("U has size: $(size(U))")
    println("D has size: $(size(D))")
    println("V has size: $(size(V))")
    
    ap_t = U * diagm(sqrt(D)) * randn(Nz, 1)
    
    # Draw a sequence of shocks, η+
    η_all_plus = sqrt(QQ)*randn(Ne,Nt0+Nt)
    
    # Set n_ant_shocks shocks to 0 in pre-ZLB time periods
    if n_ant_shocks > 0
        # get the indices of the anticipated shocks in the m.exogenous_shocks field
        ant1_ind = m.exogenous_shocks[:rm_shl1]
        antn_ind = m.exogenous_shocks[symbol("rm_shl$(n_ant_shocks)")]

        # set shocks to 0
        η_all_plus[ant1_ind:antn_ind, 1:zlb_start-1] = 0
    end
    
    # Produce "fake" states and observables (a+ and y+) by
    # iterating the state-space system forward
    for t = 1:Nt0+Nt
        ap_t             = TTT * ap_t + RRR * η_all_plus[:,t]
        α_all_plus[:,t]  = ap_t
        YY_all_plus[:,t] = ZZ*ap_t + DD
    end
    
    # Replace fake data with NaNs wherever actual data has NaNs
    YY_all = [YY0 YY];
    YY_all_plus[isnan(YY_all)] = NaN
    
    # Compute y* = y - y+ - D
    YY_star = YY_all - YY_all_plus
    
    ## Run the kalman filter

    A0, P0, pred, vpred, TTT, RRR, CCC, QQ, ZZ, DD = if n_ant_shocks > 0
        
        R2, R3, R1 = kalman_filter_2part(m, YY', allout=true)
        
        # unpack the results to pass to kalman_smoother

        TTT = R3[:TTT]
        RRR = R3[:RRR]
        CCC = R3[:CCC] 
        filtered_states = [R2[:filt] R3[:filt]]
        pred            = hcat(R2[:pred], R3[:pred])
        vpred           = cat(3, R2[:vpred], R3[:vpred])
        zend            = R3[:zend]    # final state vector is in R3
        PO_small        = R1[:P0]
        A0_small        = R1[:A0]

        r_tl1      = m.keys[:rm_tl1]
        r_tlx      = r_tl1 + 1

        # get expanded version of time-0 state vector
        A0              =  zeros(Nz)
        A0[1:r_tlx-2]   = A0[1:r_tlx-2] 
        A0[r_tl1+n_ant_shocks:end]   = A0_small[r_tlx-1:end]
                
        # get expanded version of P0 matrix
        P0              = zeros(Nz, Nz)

        P0[1:r_tlx-2, 1:r_tlx-2]                  = P0_small[1:r_tlx-2, 1:r_tlx-2]
        P0[1:r_tlx-2, (r_tl1+n_ant_shocks):end]   = P0_small[1:r_tlx-2, r_tlx-1:end]
        P0[r_tl1+n_ant_shocks:end, 1:r_tlx-2]     = P0_small[r_tlx-1:end, 1:r_tlx-2]
        P0[r_tl1+n_ant_shocks:end, r_tl1+nant:end]= P0_small[r_tlx-1:end, r_tlx-1:end]

        A0, P0, pred, vpred, TTT, RRR, CCC, QQ, ZZ, DD
    else
        myvar = zeros(Ny+Nz,Ny+Nz)
        myvar[1:Nz,1:Nz] = RRR*QQ*RRR'
        kal = kalman_filter(YY', 0, CCC, TTT, DD, ZZ, myvar, A0, P0, allout=true)

        A0, P0, kal[:pred], kal[:vpred], TTT, RRR, CCC, QQ, ZZ, DD
    end

    ##### Step 2: Kalman smooth over everything
    α_hat_star,η_hat_star = kalman_smoother(A0,P0,YY',pred,vpred,
                                            TTT,RRR,CCC,QQ,ZZ,DD,
                                            n_ant_shocks,n_ant_lags)
    
    ## Compute draw (states and shocks)
    α_til = α_plus + α_hat_star
    η_til = η_plus + η_hat_star
    
    return α_til, η_til
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
