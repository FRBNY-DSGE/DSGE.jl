# KALCVF The Kalman filter
#
# State space model is defined as follows:
#   z(t+1) = a+F*z(t)+η(t)     (state or transition equation)
#     y(t) = b+H*z(t)+ε(t)     (observation or measurement equation)
#
# [logl, <pred, vpred, <filt, vfilt>>] = kalcvf(data, lead, a, F, b, H, var, <z0, vz0>)
# computes the one-step prediction and the filtered estimate, as well as their covariance matrices.
# The function uses forward recursions, and you can also use it to obtain k-step estimates.
#
# The inputs to the KALCVF function are as follows:
#   data is a [Ny x T] matrix containing data (y(1), ... , y(T)).
#   lead is the number of steps to forecast after the end of the data.
#      a is an [Nz x 1] vector for a time-invariant input vector in the transition equation.
#      F is an [Nz x Nz] matrix for a time-invariant transition matrix in the transition equation.
#      b is an [Ny x 1] vector for a time-invariant input vector in the measurement equation.
#      H is an [Ny x Nz] matrix for a time-invariant measurement matrix in the measurement equation.
#    var is an [Ny + Nz] x [Ny + Nz] matrix for a time-invariant variance matrix for
#           the error in the transition equation and the error in the measurement equation,
#           that is, [η(t)', ε(t)']'.
#     z0 is an optional [Nz x 1] initial state vector.
#    vz0 is an optional [Nz x Nz] covariance matrix of an initial state vector.
#    Ny0 is an optional scalar indicating the number of periods of presample
#       (i.e. the number of periods which we don't add to the likelihood)
# allout is an optional keyword argument indicating whether we want optional output variables returned as well
#
# The KALCVF function returns the following output:
#   logl is a value of the average log likelihood function of the SSM
#           under assumption that observation noise ε(t) is normally distributed
#   pred is an optional [Nz x (T+lead)] matrix containing one-step predicted state vectors.
#  vpred is an optional [Nz x Nz x(T+lead)] matrix containing mean square errors of predicted state vectors.
#   filt is an optional [Nz x T] matrix containing filtered state vectors.
#  vfilt is an optional [Nz x Nz x T] matrix containing mean square errors of filtered state vectors.
#
#
# This is a M-file for MATLAB.
# Copyright 2002-2003 Federal Reserve Bank of Atlanta
# Revision: 1.2    Date: 2003/03/19 19:16:17
# Iskander Karibzhanov 5-28-02.
# Master of Science in Computational Finance
# Georgia Institute of Technology
#==========================================================================#
# Revision history:
#
#  03/19/2003  -  algorithm and interface were adapted from SAS/IML KALCVF subroutine for use in MATLAB M file
#
#==========================================================================#
function kalcvf2NaN{S<:FloatingPoint}(data::Matrix{S}, lead::Int64, a::Matrix{S}, F::Matrix{S}, b::Matrix{S}, H::Matrix{S}, var::Matrix{S}, z0::Matrix{S}, vz0::Matrix{S}, Ny0::Int = 0; allout::Bool = false)
    T = size(data, 2)
    Nz = size(a, 1)
    Ny = size(b, 1)

    z = z0
    P = vz0

    # Check input matrix dimensions
    @assert size(data, 1) == Ny
    @assert size(a, 2) == 1
    @assert size(F) == (Nz, Nz)
    @assert size(b, 2) == 1
    @assert size(H) == (Ny, Nz)
    @assert size(var) == (Ny + Nz, Ny + Nz)
    @assert size(z) == (Nz, 1)
    @assert size(P) == (Nz, Nz)

    # V(t) and R(t) are variances of η(t) and ε(t), respectively, and G(t) is a covariance of η(t) and ε(t)
    # In dsgelh :
    # --- V is same as QQ
    # --- R is same as EE
    # --- G is same as VV = QQ*MM
    V = var[1:Nz, 1:Nz]
    R = var[(Nz+1):end, (Nz+1):end]
    G = var[1:Nz, (Nz+1):end]

    if allout
        pred = zeros(Nz, T)
        vpred = zeros(Nz, Nz, T)
        
        yprederror = NaN*zeros(Ny, T)
        ystdprederror = NaN*zeros(Ny, T)
        
        filt = zeros(Nz, T)
        vfilt = zeros(Nz, Nz, T)
    end
    
    L = 0.0
    
    for t = 1:T
        # If an element of the vector y(t) is missing (NaN) for the observation t, the
        #   corresponding row is ditched from the measurement equation.
        nonmissing = !isnan(data[:, t])
        data_t = data[nonmissing, t]       # data_t = Y_T = [y1, y2, ..., yT] is matrix of observable data time-series
        H_t = H[nonmissing, :]             # H_t = DD is matrix mapping states to observables
        G_t = G[:, nonmissing]             # G_t = Cov(η_t, ε_t)
        R_t = R[nonmissing, nonmissing]    # R_t = Var(ε_t)
        Ny_t = length(data_t)              # Ny_t = T is length of time
        b_t = b[nonmissing, :]             # b_t = DD
        
        
        ## forecasting
        z = a + F*z                        # z_{t|t-1} = a + F(Θ)*z_{t-1|t-1}
        P = F*P*F' + V                     # P_{t|t-1} = F(Θ)*P_{t-1|t-1}*F(Θ)' + F(Θ)*Var(η_t)*F(Θ)'
        dy = data_t - H_t*z - b_t          # dy = y_t - H(Θ)*z_{t|t-1} - DD is prediction error or innovation
        HG = H_t*G_t                       # HG is ZZ*Cov(η_t, ε_t)
        D = H_t*P*H_t' + HG + HG' + R_t    # D = ZZ*P_{t+t-1}*ZZ' + HG + HG' + R_t 
        D = (D+D')/2

        if allout
            pred[:, t] = z
            vpred[:, :, t] = P
            yprederror[nonmissing, t] = dy
            ystdprederror[nonmissing, t] = dy./sqrt(diag(D))
        end

        ddy = D\dy 
        # We evaluate the log likelihood function by adding values of L at every iteration
        #   step (for each t = 1,2,...T)
        if t > Ny0
            L += -log(det(D))/2 - first(dy'*ddy/2) - Ny_t*log(2*pi)/2
        end
        
        ## updating
        PHG = P*H_t' + G_t
        z = z + PHG*ddy                    # z_{t|t} = z_{t|t-1} + P_{t|t-1}*H(Θ)' + ...
        P = P - PHG/D*PHG'                 # P_{t|t} = P_{t|t-1} - PHG*(1/D)*PHG

        if allout
            PH = P*H_t'
            filt[:, t] = z
            vfilt[:, :, t] = P
        end
    end
    
    zend = z
    Pend = P

    if allout && lead > 1
        for t = (T+2):(T+lead)
            z = F*z + a
            P = F*P*F' + V
            pred[:, t] = z
            vpred[:, :, t] = P
        end
    end

    if allout
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))
        return L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt
    else
        return L, zend, Pend
    end

end





# The initial state vector and its covariance matrix of the time invariant Kalman filters
# are computed under the stationarity condition:
#        z0 = (I-F)\a
#       vz0 = (I-kron(F,F))\(V(:),Nz,Nz)
# where F and V are the time invariant transition matrix and the covariance matrix of transition equation noise,
# and vec(V) is an [Nz^2 x 1] column vector that is constructed by the stacking Nz columns of matrix V.
# Note that all eigenvalues of the matrix F are inside the unit circle when the SSM is stationary.
# When the preceding formula cannot be applied, the initial state vector estimate is set to a
# and its covariance matrix is given by 1E6I. Optionally, you can specify initial values.
function kalcvf2NaN{S<:FloatingPoint}(data::Matrix{S}, lead::Int64, a::Matrix{S}, F::Matrix{S}, b::Matrix{S}, H::Matrix{S}, var::Matrix{S}, Ny0::Int = 0; allout::Bool = false)
    Nz = size(a, 1)
    V = var[1:Nz, 1:Nz]

    e, _ = eig(F)
    if countnz(e*e' - eye(Nz)) == Nz^2
        z0 = (eye(Nz) - F)\a
        vz0 = reshape((eye(Nz^2)-kron(F,F))\V, Nz, Nz)
    else
        z0 = a
        vz0 = eye(Nz)*1e6
    end

    return kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0, Ny0; allout=allout)
end

#=

# KALSMTH_K93.M

# This is a Kalman Smoothing program based on S.J. Koopman's "Disturbance
# Smoother for State Space Models" (Biometrika, 1993), as specified in
# Durbin and Koopman's "A Simple and Efficient Simulation Smoother for
# State Space Time Series Analysis" (Biometrika, 2002). The algorithm has been
# simplified for the case in which there is no measurement error, and the
# model matrices do not vary with time.

# Unlike other Kalman Smoothing programs, there is no need to invert
# singular matrices using the Moore-Penrose pseudoinverse (pinv), which
# should lead to efficiency gains and fewer inversion problems. Also, the
# states vector and the corresponding matrices do not need to be augmented
# to include the shock innovations. Instead they are saved automatically
# in the eta_hat matrix.

# Nz will stand for the number of states, Ny for the number of observables,
# Ne for the number of shocks, and Nt for the number of periods of data.

# The state space is assumed to take the form:
# y(t) = Z*alpha(t) + b
# alpha(t+1) = T*alpha(t) + R*eta(t+1)

# INPUTS:

# A0, the (Nz x 1) initial (time 0) states vector.
# P0, the (Nz x Nz) initial (time 0) state covariance matrix.
# y, the (Ny x Nt) matrix of observable data.
# pred, the (Nz x Nt) matrix of one-step-ahead predicted states (from the Kalman Filter).
# vpred, the (Nz x Nz x Nt) matrix of one-step-ahead predicted covariance matrices.
# T, the (Nz x Nz) transition matrix.
# R, the (Nz x Ne) matrix translating shocks to states.
# Q, the (Ne x Ne) covariance matrix for the shocks.
# Z, the (Ny x Nz) measurement matrix.
# b, the (Ny x 1) constant vector in the measurement equation.

# nant, an optional scalar for the zero bound specification indicating the
#       number of periods ahead the interest rate is fixed.
# antlags, an optional scalar for the zero bound specification indicating
#       the number of periods for which interest rate expectations have
#       been fixed
# Ny0, an optional scalar indicating the number of periods of presample
#       (i.e. the number of periods for which smoothed states are not required).

# OUTPUTS:

# alpha_hat, the (Nz x Nt) matrix of smoothed states.
# eta_hat, the optional (Ne x Nt) matrix of smoothed shocks.

# If Ny0 is nonzero, the alpha_hat and eta_hat matrices will be shorter by
# that number of columns (taken from the beginning).

# Dan Greenwald, 7/7/2010.

function kalsmth_k93(A0, P0, y, pred::Matrix{S}, vpred::Array{S, 3}, T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S}, b::Matrix{S}, nant::Int, antlags::Int, peachcount, psize, Ny0::Int = 0)

    Ne = size(R, 2)
    Ny = size(y, 1)
    Nt = size(y, 2)
    Nz = size(T, 1)

    alpha_hat = zeros(Nz, Nt)

    [r, eta_hat] = distsmth_k93(y, pred, vpred, T, R, Q, Z, b, peachcount, psize, nant, antlags)
    
    ah_t = A0 + P0*r[:, 1]
    alpha_hat[:, 1] = ah_t

    for t = 2:Nt

        # This section relates to the zero bound framework, in which no
        # anticipated shocks are supposed to occur before the model switch.
        # In these periods, this is accomplished by setting the relevant
        # rows and columns of the Q matrix to zero. In other periods, or in
        # specifications with zero bound off (and hence with nant = 0), the
        # normal Q matrix can be used.

        if nant != 0
            # The first part of the conditional below pertains to the periods in which zerobound is off.
            # To specify this period, we must account for (peachcount*psize) since peachdata is augmented to y.
            # JC 11/30/10
            if nant > 0 && t < Nt - antlags - (peachcount*psize)
                Q_t = zeros(Ne, Ne)
                Q_t[1:(Ne-nant), 1:(Ne-nant)] = Q[1:(Ne-nant), 1:(Ne-nant)]
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

# nant, an optional scalar for the zero bound specification indicating the
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
function distsmth_k93{S<:FloatingPoint}(y::Matrix{S}, pred::Matrix{S}, vpred::Matrix{S}, T::Matrix{S}, R::Matrix{S}, Q::Matrix{S}, Z::Matrix{S}, b::Matrix{S}, peachcount, psize, nant::Int = 0, antlags::Int = 0)
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
        # specifications with zero bound off (and hence with nant = 0), the
        # normal Q matrix can be used.
        if nant != 0
            
            # The first part of the conditional below pertains to the periods in which zerobound is off.
            # To specify this period, we must account for (peachcount*psize) since peachdata is augmented to y.
            # JC 11/30/10
            if nant > 0 && t < Nt - antlags - peachcount*psize
                Q_t = zeros(Ne, Ne)
                Q_t[1:Ne-nant, 1:Ne-nant] = Q[1:Ne-nant, 1:Ne-nant]
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
=#
