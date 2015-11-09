
"""
`kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0, Ny0; allout=false)`
`kalcvf2NaN(data, lead, a, F, b, H, var, Ny0=0; allout=false)`

Inputs
------

- `data` is a [Ny x T] matrix containing data (y(1), ... , y(T)).
- `lead` is the number of steps to forecast after the end of the data.
- `a` is an [Nz x 1] vector for a time-invariant input vector in the transition equation.
- `F` is an [Nz x Nz] matrix for a time-invariant transition matrix in the transition
  equation.
- `b` is an [Ny x 1] vector for a time-invariant input vector in the measurement equation.
- `H` is an [Ny x Nz] matrix for a time-invariant measurement matrix in the measurement
  equation.
- `var` is an [Ny + Nz] x [Ny + Nz] matrix for a time-invariant variance matrix for the
  error in the transition equation and the error in the measurement equation, that is,
  [η(t)', ϵ(t)']'.
- `z0` is an optional [Nz x 1] initial state vector.
- `vz0` is an optional [Nz x Nz] covariance matrix of an initial state vector.
- `Ny0` is an optional scalar indicating the number of periods of presample (i.e. the number
  of periods which we don't add to the likelihood)
- `allout` is an optional keyword argument indicating whether we want optional output
  variables returned as well 


Outputs
-------

- `logl` is a value of the average log likelihood function of the SSM under assumption that
  observation noise ϵ(t) is normally distributed
- `pred` is a [Nz x (T+lead)] matrix containing one-step predicted state vectors.
- `vpred` is a [Nz x Nz x(T+lead)] matrix containing mean square errors of predicted
  state vectors.
- `filt` is an optional [Nz x T] matrix containing filtered state vectors.
- `vfilt` is an optional [Nz x Nz x T] matrix containing mean square errors of filtered state
  vectors.


Notes
-----
State space model is defined as follows:
```
z(t+1) = a+F*z(t)+η(t)     (state or transition equation)
y(t) = b+H*z(t)+ϵ(t)       (observation or measurement equation)
```

When z0 and Vz0 are omitted, the initial state vector and its covariance matrix of the time
invariant Kalman filters are computed under the stationarity condition:
```
z0 = (I-F)\a
vz0 = (I-kron(F,F))\(V(:),Nz,Nz)
```
where F and V are the time invariant transition matrix and the covariance matrix of
transition equation noise, and vec(V) is an [Nz^2 x 1] column vector that is constructed by
the stacking Nz columns of matrix V.  Note that all eigenvalues of the matrix F are inside
the unit circle when the SSM is stationary.  When the preceding formula cannot be applied,
the initial state vector estimate is set to a and its covariance matrix is given by 1E6I.
Optionally, you can specify initial values.

Attribution
-----------
Adapted from `KALCVF`, Iskander Karibzhanov, Federal Reserve Bank of Atlanta, 2003-03-19.
"""
function kalcvf2NaN{S<:AbstractFloat}(data::Matrix{S},
                                      lead::Int64,
                                      a::Matrix{S},
                                      F::Matrix{S},
                                      b::Matrix{S},
                                      H::Matrix{S},
                                      var::Matrix{S},
                                      z0::Matrix{S},
                                      vz0::Matrix{S},
                                      Ny0::Int = 0;
                                      allout::Bool = false)
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

    # V(t) and R(t) are variances of η(t) and ϵ(t), respectively, and G(t) is a covariance
    # of η(t) and ϵ(t)
    # In dsgelh :
    # --- V is same as QQ
    # --- R is same as EE
    # --- G is same as VV = QQ*MM
    V = var[1:Nz, 1:Nz]
    R = var[(Nz+1):end, (Nz+1):end]
    G = var[1:Nz, (Nz+1):end]

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
        H_t = H[nonmissing, :]             # H_t = DD is matrix mapping states to observables
        G_t = G[:, nonmissing]             # G_t = Cov(η_t, ϵ_t)
        R_t = R[nonmissing, nonmissing]    # R_t = Var(ϵ_t)
        Ny_t = length(data_t)              # Ny_t = T is length of time
        b_t = b[nonmissing, :]             # b_t = DD
        
        
        ## forecasting
        z = a + F*z                        # z_{t|t-1} = a + F(Θ)*z_{t-1|t-1}
        P = F*P*F' + V                     # P_{t|t-1} = F(Θ)*P_{t-1|t-1}*F(Θ)' + F(Θ)*Var(η_t)*F(Θ)'
        dy = data_t - H_t*z - b_t          # dy = y_t - H(Θ)*z_{t|t-1} - DD is prediction error or innovation
        HG = H_t*G_t                       # HG is ZZ*Cov(η_t, ϵ_t)
        D = H_t*P*H_t' + HG + HG' + R_t    # D = ZZ*P_{t+t-1}*ZZ' + HG + HG' + R_t 
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
        if t > Ny0
            L += -log(det(D))/2 - first(dy'*ddy/2) - Ny_t*log(2*pi)/2
        end
        
        ## updating
        PHG = P*H_t' + G_t
        z = z + PHG*ddy                    # z_{t|t} = z_{t|t-1} + P_{t|t-1}*H(Θ)' + ...
        P = P - PHG/D*PHG'                 # P_{t|t} = P_{t|t-1} - PHG*(1/D)*PHG

        if allout
            PH = P*H_t'
            filt[:, t]     = z
            vfilt[:, :, t] = P
        end
    end
    
    zend = z
    Pend = P

    if allout && lead > 1
        for t = (T+2):(T+lead)
            z = F*z + a
            P = F*P*F' + V
            pred[:, t]     = z
            vpred[:, :, t] = P
        end
    end

    if allout
        rmse = sqrt(mean((yprederror.^2)', 1))
        rmsd = sqrt(mean((ystdprederror.^2)', 1))
        return Kalman(L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt)
    else
        return Kalman(L, zend, Pend)
    end

end

function kalcvf2NaN{S<:AbstractFloat}(data::Matrix{S},
                                      lead::Int64,
                                      a::Matrix{S},
                                      F::Matrix{S},
                                      b::Matrix{S},
                                      H::Matrix{S},
                                      var::Matrix{S},
                                      Ny0::Int = 0;
                                      allout::Bool = false)
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

immutable Kalman{S<:AbstractFloat}
    L::S
    zend::Matrix{S}
    Pend::Matrix{S}
    pred::Matrix{S}
    vpred::Array{S,3}
    yprederror::Matrix{S}
    ystdprederror::Matrix{S}
    rmse::Matrix{S}
    rmsd::Matrix{S}
    filt::Matrix{S}
    vfilt::Array{S,3}
end
function Kalman{S<:AbstractFloat}(L::S,
                                  zend::Matrix{S},
                                  Pend::Matrix{S},
                                  pred::Matrix{S}          = Matrix{S}(),
                                  vpred::Array{S,3}        = Array{S}(0,0,0),
                                  yprederror::Matrix{S}    = Matrix{S}(),
                                  ystdprederror::Matrix{S} = Matrix{S}(),
                                  rmse::Matrix{S}          = Matrix{S}(),
                                  rmsd::Matrix{S}          = Matrix{S}(),
                                  filt::Matrix{S}          = Matrix{S}(),
                                  vfilt::Array{S,3}        = Array{S}(0,0,0))
    return Kalman{S}(L,zend,Pend,pred,vpred,yprederror,ystdprederror,rmse,rmsd,filt,vfilt)
end
function Base.getindex(K::Kalman, d::Symbol)
    if d in (:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse, :rmsd,
             :filt, :vfilt)
        return getfield(K, d)
    else
        throw(KeyError(d))
    end
end
