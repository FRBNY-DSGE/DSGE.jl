module Kalman

export kalcvf2NaN

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
function kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0)
    T = size(data,2)
    Nz = size(a,1)
    Ny = size(b,1)

    z = z0
    P = vz0

    # Check input matrix dimensions
    if size(data, 1) != Ny
        error("data and b must have the same number of rows")
    end
    if size(a, 2) != 1
        error("a must be column vector")
    end
    if size(F) != (Nz, Nz)
        error("F must be square")
    end
    if size(b, 2) != 1
        error("b must be column vector")
    end
    if size(H) != (Ny, Nz)
        error("H must be Ny by Nz matrix")
    end
    if size(var) != (Ny+Nz, Ny+Nz)
        error("var must be (Ny+Nz) by (Ny+Nz) matrix")
    end
    if size(z) != (Nz, 1)
        error("z0 must be column vector of length Nz")
    end
    if size(P) != (Nz, Nz)
        error("vz0 must be Nz by Nz matrix")
    end

    # V(t) and R(t) are variances of η(t) and ε(t), respectively, and G(t) is a covariance of η(t) and ε(t)
    # In dsgelh :
    # --- V is same as QQ
    # --- R is same as EE
    # --- G is same as VV = QQ*MM
    V = var[1:Nz, 1:Nz]
    R = var[Nz+1:end, Nz+1:end]
    G = var[1:Nz, Nz+1:end]
    
    pred = zeros(Nz, T)
    vpred = zeros(Nz, Nz, T)
    
    yprederror = NaN*zeros(Ny, T)
    ystdprederror = NaN*zeros(Ny, T)

    filt = zeros(Nz, T)
    vfilt = zeros(Nz, Nz, T)
    
    L = 0
    
    for t = 1:T
        # If an element of the vector y(t) is missing (NaN) for the observation t, the
        #   corresponding row is ditched from the measurement equation.
        notis_nan = !isnan(data[:, t])
        data_t = data[notis_nan, t]        # data_t = Y_T = [y1, y2, ..., yT] is matrix of observable data time-series
        H_t = H[notis_nan, :]              # H_t = DD is matrix mapping states to observables
        G_t = G[:, notis_nan]              # G_t = Cov(η_t, ε_t)
        R_t = R[notis_nan, notis_nan]      # R_t = Var(ε_t)
        Ny_t = length(data_t)              # Ny_t = T is length of time
        b_t = b[notis_nan, :]              # b_t = DD
        
        
        ## forecasting
        z = a + F*z                        # z_{t|t-1} = a + F(Θ)*z_{t-1|t-1}
        P = F*P*F' + V                     # P_{t|t-1} = F(Θ)*P_{t-1|t-1}*F(Θ)' + F(Θ)*Var(η_t)*F(Θ)'
        dy = data_t - H_t*z - b_t          # dy = y_t - H(Θ)*z_{t|t-1} - DD is prediction error or innovation
        HG = H_t*G_t                       # HG is ZZ*Cov(η_t, ε_t)
        D = H_t*P*H_t' + HG + HG' + R_t    # D = ZZ*P_{t+t-1}*ZZ' + HG + HG' + R_t 
        D = .5*(D+D')

        pred[:, t] = z
        vpred[:, :, t] = P
        yprederror[notis_nan, t] = dy
        ystdprederror[notis_nan, t] = dy./sqrt(diag(D))

        ddy = D\dy 
        # We evaluate the log likelihood function by adding values of L at every iteration
        #   step (for each t = 1,2,...T)
        L += -log(det(D))/2 - dy'*ddy/2 - Ny_t*log(2*pi)/2

        
        ## updating
        PHG = P*H_t' + G_t
        z = z + PHG*ddy                    # z_{t|t} = z_{t|t-1} + P_{t|t-1}*H(Θ)' + ...
        P = P - PHG/D*PHG'                 # P_{t|t} = P_{t|t-1} - PHG*(1/D)*PHG

        PH = P*H_t'
        filt[:, t] = z
        vfilt[:, :, t] = P
    end
    
    zend = z
    Pend = P

    if lead > 1
        for t = T+2:T+lead
            z = F*z + a
            P = F*P*F' + V
            pred[:, t] = z
            vpred[:, :, t] = P
        end
    end

    rmse = sqrt(mean((yprederror.^2)', 1))
    rmsd = sqrt(mean((ystdprederror.^2)', 1))

    return L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt
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
function kalcvf2NaN(data, lead, a, F, b, H, var)
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

    return kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0)
end

end # module
