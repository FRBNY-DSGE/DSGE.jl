#=
This code is loosely based on
Strid and Walentin (2008) "Block Kalman filtering for large-scale DSGE models"
=#

mutable struct BlockKalmanFilter{S<:AbstractFloat}
    A1::Matrix{S}
    A2::Matrix{S}
    A3::Matrix{S}
    B1::Matrix{S}
    B2::Matrix{S}
    B3::Matrix{S}
    Cblock::Matrix{S}
    Rup::Matrix{S}
    Rlo::Matrix{S}
    C::Vector{S}
    Q::Matrix{S}
    Z1::Matrix{S}
    Z2::Matrix{S}
    Z3::Matrix{S}
    Z4::Matrix{S}
    D::Vector{S}
    E::Matrix{S}
    s1_t::Vector{S} # s1_{t|t-1} or s1_{t|t}, and so on
    s2_t::Vector{S}
    s3_t::Vector{S}
    s4_t::Vector{S}
    P11_t::Matrix{S} # P11_{t|t-1} or P11_{t|t}
    P12_t::Matrix{S}
    P13_t::Matrix{S}
    P14_t::Matrix{S}
    P22_t::Matrix{S}
    P23_t::Matrix{S}
    P24_t::Matrix{S}
    P33_t::Matrix{S}
    P34_t::Matrix{S}
    P44_t::Matrix{S}
    M::Matrix{S}
    Mtild::Matrix{S}
    block_dims::Vector{Int64}
    block_num::Int64
    true_block::Bool
    loglh_t::S     # P(y_t | y_{1:t})
end

"""
```
KalmanFilter(T, R, C, Q, Z, D, E, [s_0, P_0])
```

Outer constructor for the `KalmanFilter` type.
"""
function BlockKalmanFilter(Ttild::Matrix{S}, Rtild::Matrix{S}, Ctild::Vector{S}, Qtild::Matrix{S},
                      Ztild::Matrix{S}, Dtild::Vector{S}, Etild::Matrix{S}, M::Matrix{S}, Mtild::Matrix{S},
                           block_dims::Vector{Int64}, s_0tild::Vector{S} = Vector{S}(0),
                           P_0tild::Matrix{S} = Matrix{S}(0, 0); block_num::Int64 = 2, true_block::Bool = false) where {S<:AbstractFloat}
    if isempty(s_0tild) || isempty(P_0tild)
        s_0tild, P_0tild = init_stationary_states(Ttild, Rtild, Ctild, Qtild)
    end
    T = M' * Ttild * M
    R = M' * Rtild * Mtild
    C = M * Ctild
    Q = Mtild' * Qtild * Mtild
    Z = Ztild * M
    # D and E are left as is
    s_0 = M * s_0tild
    P_0 = M' * P_0tild * M
    if block_dims[2] == 0 & block_dims[3] == 0
        # 2-block case
        dim1 = block_dims[1]; dim2 = block_dims[4]
        A1 = T[1:dim1, 1:dim1]
        B1 = R[dim1 + 1:end, 1:dim1]
        Cblock = T[dim1 + 1:end, dim1 + 1:end]
        Z1 = Matrix{S}(0, 0)
        Z2 = Z[:, dim1 + 1:end]
        Rup = Matrix{S}(0, 0)
        Rlo = Matrix{S}(0, 0)

        # Check if dim1 is equal to number of exogenous shocks; if not, need to track more matrices
        if dim1 != size(R, 2)
            true_block = false
            Rup = R[1:dim1, dim1 + 1:end]
            Rlo = R[dim1 + 1:end, dim1 + 1:end]
        end
        # Check Z is in true block form; if not (e.g. unit root), need to track more matrices
        if Z[:, 1:dim1] != zeros(Z[:, 1:dim1])
            true_block = false
            Z1 = Z[:, 1:dim1]
        end
        s1_0 = s_0[1:dim1]
        s2_0 = s_0[dim1 + 1:end]
        P11_0 = P_0[1:dim1, 1:dim1]
        P12_0 = P_0[1:dim1, dim1 + 1:end]
        P22_0 = P_0[dim1 + 1:end, dim1 + 1:end]

        return BlockKalmanFilter(A1, Matrix{S}(0, 0), Matrix{S}(0, 0), B1, Matrix{S}(0, 0), Matrix{S}(0, 0),
                          Cblock, Rup, Rlo, C, Q, Z1, Z2, Matrix{S}(0, 0),
                          Matrix{S}(0, 0), D, E, s1_0, s2_0, Vector{S}(0), Vector{S}(0),
                          P11_0, P12_0, Matrix{S}(0, 0), Matrix{S}(0, 0), P22_0,
                          Matrix{S}(0, 0), Matrix{S}(0, 0), Matrix{S}(0, 0), Matrix{S}(0, 0),
                          Matrix{S}(0, 0), M, Mtild, block_dims, block_num, true_block, NaN)
    end
end

"""
```
init_stationary_states(T, R, C, Q)
```

Compute the initial state `s_0` and state covariance matrix `P_0` under the
stationarity condition:

```
s_0  = (I - T)\C
P_0 = reshape(I - kron(T, T))\vec(R*Q*R'), Ns, Ns)
```

where:

- `kron(T, T)` is a matrix of dimension `Ns^2` x `Ns^2`, the Kronecker
  product of `T`
- `vec(R*Q*R')` is the `Ns^2` x 1 column vector constructed by stacking the
  `Ns` columns of `R*Q*R'`

All eigenvalues of `T` are inside the unit circle when the state space model
is stationary. When the preceding formula cannot be applied, the initial state
vector estimate is set to `C` and its covariance matrix is given by `1e6 * I`.
"""
function init_stationary_states(T::Matrix{S}, R::Matrix{S}, C::Vector{S},
                                Q::Matrix{S}) where {S<:AbstractFloat}
    e, _ = eig(T)
    if all(abs.(e) .< 1)
        s_0 = (UniformScaling(1) - T)\C
        P_0 = solve_discrete_lyapunov(T, R*Q*R')
    else
        Ns = size(T, 1)
        s_0 = C
        P_0 = 1e6 * eye(Ns)
    end
    return s_0, P_0
end

"""
```
kalman_filter(y, T, R, C, Q, Z, D, E, s_0 = Vector(), P_0 = Matrix();
    outputs = [:loglh, :pred, :filt], Nt0 = 0)

kalman_filter(regime_indices, y, Ts, Rs, Cs, Qs, Zs, Ds, Es,
    s_0 = Vector(), P_0 = Matrix(); outputs = [:loglh, :pred, :filt],
    Nt0 = 0)
```

This function implements the Kalman filter for the following state-space model:

```
s_{t+1} = C + T*s_t + R*ϵ_t    (transition equation)
y_t     = D + Z*s_t + u_t      (measurement equation)

ϵ_t ∼ N(0, Q)
u_t ∼ N(0, E)
Cov(ϵ_t, u_t) = 0
```

### Inputs

- `y`: `Ny` x `Nt` matrix containing data `y_1, ... , y_T`
- `s_0`: optional `Ns` x 1 initial state vector
- `P_0`: optional `Ns` x `Ns` initial state covariance matrix

**Method 1 only:**

- `T`: `Ns` x `Ns` state transition matrix
- `R`: `Ns` x `Ne` matrix in the transition equation mapping shocks to states
- `C`: `Ns` x 1 constant vector in the transition equation
- `Q`: `Ne` x `Ne` matrix of shock covariances
- `Z`: `Ny` x `Ns` matrix in the measurement equation mapping states to
  observables
- `D`: `Ny` x 1 constant vector in the measurement equation
- `E`: `Ny` x `Ny` matrix of measurement error covariances

**Method 2 only:**

- `regime_indices`: `Vector{Range{Int}}` of length `n_regimes`, where
  `regime_indices[i]` indicates the time periods `t` in regime `i`
- `Ts`: `Vector{Matrix{S}}` of `T` matrices for each regime
- `Rs`
- `Cs`
- `Qs`
- `Zs`
- `Ds`
- `Es`

where:

- `Nt`: number of time periods for which we have data
- `Ns`: number of states
- `Ne`: number of shocks
- `Ny`: number of observables

### Keyword Arguments

- `outputs`: some subset of `[:loglh, :pred, :filt]` specifying which outputs to
  compute and return. There will always be the same number of return values,
  but, for example, `s_pred` and `P_pred` will be returned as empty arrays if
  `:pred` is not in `outputs`
- `Nt0`: number of presample periods to omit from all return values

### Outputs

- `loglh`: length `Nt` vector of conditional log-likelihoods P(y_t | y_{1:t-1})
- `s_pred`: `Ns` x `Nt` matrix of one-step predicted state vectors s_{t|t-1}
- `P_pred`: `Ns` x `Ns` x `Nt` array of mean squared errors P_{t|t-1} of
  predicted state vectors
- `s_filt`: `Ns` x `Nt` matrix of filtered state vectors s_{t|t}
- `P_filt`: `Ns` x `Ns` x `Nt` matrix containing mean squared errors P_{t|t} of
  filtered state vectors
- `s_0`: `Ns` x 1 initial state vector. This may have been reassigned to the
  last presample state vector if `Nt0 > 0`
- `P_0`: `Ns` x `Ns` initial state covariance matrix. This may have been
  reassigned to the last presample state covariance if `Nt0 > 0`
- `s_T`: `Ns` x 1 final filtered state `s_{T|T}`
- `P_T`: `Ns` x `Ns` final filtered state covariance matrix `P_{T|T}`

### Notes

When `s_0` and `P_0` are omitted, they are computed using
`init_stationary_states`.
"""
function block_kalman_filter(y::Matrix{Float64}, Ttild::Matrix{Float64}, Rtild::Matrix{Float64},
                             Ctild::Vector{Float64}, Qtild::Matrix{Float64}, Ztild::Matrix{Float64},
                             Dtild::Vector{Float64}, Etild::Matrix{Float64},
                             M::Matrix{Float64}, Mtild::Matrix{Float64}, block_dims::Vector{Int64},
                             s_0tild::Vector{Float64} = Vector{Float64}(0),
                             P_0tild::Matrix{Float64} = Matrix{Float64}(0,0);
                             outputs::Vector{Symbol} = [:loglh, :pred, :filt],
                             Nt0::Int = 0, block_num::Int64 = 2, true_block::Bool = false)

    # Determine outputs
    return_loglh = :loglh in outputs
    return_pred  = :pred in outputs
    return_filt  = :filt in outputs

    # Dimensions
    Ns = size(T,1) # number of states
    Nt = size(y,2) # number of periods of data

    # Initialize Inputs and outputs, populate initial states
    s_0 = M * s_0tild
    P_0 = M' * P_0tild * M
    k = BlockKalmanFilter(Ttild, Rtild, Ctild, Qtild, Ztild, Dtild, Etild, M, Mtild,
                          block_dims, s_0tild, P_0tild; block_num = block_num, true_block = true_block)

    mynan  = convert(Float64, NaN)
    loglh  = return_loglh ? fill(mynan, Nt)         : Vector{Float64}(0)
    s_pred = return_pred  ? fill(mynan, Ns, Nt)     : Matrix{Float64}(0, 0)
    P_pred = return_pred  ? fill(mynan, Ns, Ns, Nt) : Array{Float64, 3}(0, 0, 0)
    s_filt = return_filt  ? fill(mynan, Ns, Nt)     : Matrix{Float64}(0, 0)
    P_filt = return_filt  ? fill(mynan, Ns, Ns, Nt) : Array{Float64, 3}(0, 0, 0)

    # Loop through periods t
    for t = 1:Nt
        # Forecast
        forecast!(k)
        if return_pred
            if k.block_num == 2
                s_pred[:,    t] = vcat(k.s1_t, k.s2_t)
                P_pred[:, :, t] = [k.P11_t k.P12_t; k.P12_t' k.P22_t]
            end
        end

        # Update and compute log-likelihood
        update!(k, y[:, t]; return_loglh = return_loglh)
        if return_filt
            if k.block_num == 2
                s_filt[:,    t] = vcat(k.s1_t, k.s2_t)
                P_filt[:, :, t] = [k.P11_t k.P12_t; k.P12_t' k.P22_t]
            end
        end
        if return_loglh
            loglh[t]        = k.loglh_t
        end

        # # Update s_0 and P_0 if Nt_0 > 0, presample periods
        # if t == Nt0
        #     s_0 = k.s_t
        #     P_0 = k.P_t
        # end
    end

    # Populate final states
    if k.block_num == 2
        s_T = vcat(k.s1_t, k.s2_t)
        P_T = [k.P11_t k.P12_t; k.P12_t' k.P22_t]
    end

    # # Remove presample periods
    # loglh, s_pred, P_pred, s_filt, P_filt =
    #     remove_presample!(Nt0, loglh, s_pred, P_pred, s_filt, P_filt; outputs = outputs)

    return loglh, s_pred, P_pred, s_filt, P_filt, s_0, P_0, s_T, P_T
end

# Computes the one-step-ahead states s_{t|t-1} and state covariances P_{t|t-1}
# and assign to 'k'
function forecast!(k::BlockKalmanFilter{Float64})

    if k.block_num == 2
        # Get block matrices that don't depend on true_block
        dim1 = k.block_dims[1] # dimension of exogenous AR(1) processes
        dim2 = k.block_dims[4] # dimension of endogenous states
        s1_filt = k.s1_t
        s2_filt = k.s2_t
        P11_filt = k.P11_t
        P12_filt = k.P12_t
        P22_filt = k.P22_t
        Ablock = k.A1
        Bblock = k.B1
        Cblock = k.Cblock
        Q1 = k.Q[1:dim1, 1:dim1]
        Q2 = k.Q[dim1 + 1:end, dim1 + 1:end]

        # Compute P_t, s_t
        if k.true_block
            L_t = (P12_filt * Cblock') .* (diag(Ablock) * ones(1, dim2))
            P11_t = P11_filt .* (diag(Ablock) * diag(Ablock)') + Q1
            P12_t = P11_t * Bblock' + L_t
            G = Bblock * (P12_t + L_t)
            P22_t = (G + G')/2 + Cblock * P22_filt * Cblock'

            s1_t = diag(Ablock) .* s1_filt
            s2_t = Bblock * k.s1_t + Cblock * s2_filt
        else
            # RQR' = |  Q1 + Rup * Q2 * Rup', Q1*Bblock' + Rup * Q2 * Rlo'        |
            #        |     ---           , Bblock * Q1 * Bblock' + Rlo * Q2 * Rlo'|
            L_t = (P12_filt * Cblock') .* (diag(Ablock) * ones(1, dim2))
            RupQ2 = k.Rup * Q2

            # Compute intermediate matrices to minimize additions
            P11_t = P11_filt .* (diag(Ablock) * diag(Ablock)') + Q1
            P12_t = P11_t * Bblock' + L_t
            G = Bblock * (P12_t + L_t)

            # Compute final matrices
            P22_t = (G + G')/2 + k.Rlo * Q2 * k.Rlo' + Cblock * P22_filt * Cblock'
            P12_t = P12_t + RupQ2 * k.Rlo'
            P11_t = P11_t + RupQ2 * k.Rup'

            s1_t = diag(Ablock) .* s1_filt
            s2_t = Bblock * s1_t + Cblock * s2_filt
        end

        # Save values
        k.s1_t = s1_t; k.s2_t = s2_t; k.P11_t = P11_t; k.P12_t = P12_t; k.P22_t = P22_t
    end
    return nothing
end

function update!(k::BlockKalmanFilter{Float64}, y_obs::Vector{Float64}; return_loglh::Bool = true)
    # Keep rows of measurement equation corresponding to non-NaN observables
    nonnan = .!isnan.(y_obs)
    y_obs = y_obs[nonnan]
    D = k.D[nonnan]
    E = k.E[nonnan, nonnan]
    Ny = length(y_obs)

    if k.block_num == 2
        # Get block matrices
        dim1 = k.block_dims[1] # dimension of exogenous AR(1) processes
        dim2 = k.block_dims[4] # dimension of endogenous states
        s1_pred  = k.s1_t
        s2_pred  = k.s2_t
        P11_pred = k.P11_t
        P12_pred = k.P12_t
        P22_pred = k.P22_t
        Z2 = k.Z2[nonnan, :]

        # Compute predicted y, measurement covariance matrix, error, and auxiliary matrices
        if k.true_block
            y_pred = Z2 * s2_pred + D
            V_pred = Z2 * P22_pred * Z2' + E
        else
            Z1 = k.Z1[nonnan, :]
            off_diag_mat = Z1 * P12_pred * Z2'
            y_pred = Z1 * s1_pred + Z2 * s2_pred + D
            V_pred = Z1 * P11_pred * Z1' + off_diag_mat + off_diag_mat' + Z2 * P22_pred * Z2' + E
        end
        V_pred_inv = inv(V_pred)
        dy = y_obs - y_pred # prediction error

        # Update state covariance matrix
        if k.true_block
            PZV_1 = P12_pred * Z2' * V_pred_inv
            PZV_2 = P22_pred * Z2' * V_pred_inv

            Z2_P22_pred = Z2 * P22_pred
            P11_t = P11_pred - PZV_1 * Z2 * P12_pred'
            P12_t = P12_pred - PZV_1 * Z2_P22_pred
            P22_t = P22_pred - PZV_2 * Z2_P22_pred
        else
            kalman_gain_1 = P11_pred * Z1' + P12_pred * Z2'
            kalman_gain_2 = P12_pred' * Z1' + P22_pred * Z2'
            PZV_1 = kalman_gain_1 * V_pred_inv
            PZV_2 = kalman_gain_2 * V_pred_inv

            P11_t = P11_pred - PZV_1 * kalman_gain_1'
            P12_t = P12_pred - PZV_1 * kalman_gain_2'
            P22_t = P22_pred - PZV_2 * kalman_gain_2'
        end

        # Update states
        s1_t = s1_pred + PZV_1 * dy
        s2_t = s2_pred + PZV_2 * dy

        # Save matrices
        k.s1_t = s1_t; k.s2_t = s2_t
        k.P11_t = P11_t; k.P12_t = P12_t; k.P22_t = P22_t

        if return_loglh
            k.loglh_t = -(Ny * log(2*π) + log(det(V_pred)) + dy'*V_pred_inv*dy)/2
        end
    end
end

# """
# ```
# remove_presample!(Nt0, loglh, s_pred, P_pred, s_filt, P_filt)
# ```

# Remove the first `Nt0` periods from all other input arguments and return.
# """
# function remove_presample!(Nt0::Int, loglh::Vector{S},
#                            s_pred::Matrix{S}, P_pred::Array{S, 3},
#                            s_filt::Matrix{S}, P_filt::Array{S, 3};
#                            outputs::Vector{Symbol} = [:loglh, :pred, :filt]) where {S<:AbstractFloat}
#     if Nt0 > 0
#         if :loglh in outputs
#             loglh  = loglh[(Nt0+1):end]
#         end
#         if :pred in outputs
#             s_pred = s_pred[:,    (Nt0+1):end]
#             P_pred = P_pred[:, :, (Nt0+1):end]
#         end
#         if :filt in outputs
#             s_filt = s_filt[:,    (Nt0+1):end]
#             P_filt = P_filt[:, :, (Nt0+1):end]
#         end
#     end
#     return loglh, s_pred, P_pred, s_filt, P_filt
# end
