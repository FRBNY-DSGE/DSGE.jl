#=
This code is loosely based on
Strid and Walentin (2008) "Block Kalman filtering for large-scale DSGE models"
=#

mutable struct CTBlockKalmanFilter{S<:AbstractFloat}
    T_powers::Dict{Int64, Matrix{S}}
    R_blocks::Dict{Int64, Matrix{S}}
    C::Vector{S}
    Q::Matrix{S}
    Z::Matrix{S}
    D::Vector{S}
    E::Matrix{S}
    RQRp::Dict{Tuple{Int64, Int64}, Matrix{S}} # RQR', indexed block-wise
    s_t::Vector{S}
    P_t::Dict{Tuple{Int64, Int64}, Matrix{S}}  # P_t indexed block-wise
    n_simulate_states::Int64
    loglh_t::S     # P(y_t | y_{1:t})
end

"""
```
KalmanFilter(T, R, C, Q, Z, D, E, [s_0, P_0])
```

Outer constructor for the `KalmanFilter` type.
"""
function CTBlockKalmanFilter(T::Matrix{S}, R::Matrix{S}, C::Vector{S}, Q::Matrix{S},
                             Z::Matrix{S}, D::Vector{S}, E::Matrix{S}, s_0::Vector{S} = Vector{S}(0),
                             P_0::Matrix{S} = Matrix{S}(0, 0),
                             n_simulate_states::Int64=1) where {S<:AbstractFloat}
    if isempty(s_0) || isempty(P_0)
        s_0, P_0 = init_stationary_states(T, R, C, Q)
    end

    # Create dictionaries of products of T matrices, plus create blocks of R matrices/extend C -> TT, RR, CC
    T_powers = Dict{Int64, Matrix{S}}()  # key is power of T
    R_blocks = Dict{Int64, Matrix{S}}() # key is power of T multiplying R
    R_blocks[0] = R
    T_powers[1] = T
    for i = 2:n_simulate_states
        T_powers[i] = T_powers[i - 1] * T
        R_blocks[i - 1] = T_powers[i - 1] * R
    end

    # Create RQR' dictionary
    RQRp[(1,1)] = R_blocks[0] * Q * R_blocks[0]' # initialize first row as base case
    if n_simulate_states > 1
        for i = 2:n_simulate_states
            RQRp[(1,i)] = R-blocks[0] * Q * R_blocks[i - 1]'
        end

        # Recursive construction of rest of RQR' dictionary
        for i = 2:n_simulate_states
            for j = 1:n_simulate_states
                if i <= j
                    RQRp[(i, j)] = R_blocks[i - 1] * Q * R_blocks[j - 1] + RQRp[(i-1,j-1)]
                else
                    # RQR' is symmetric
                    RQRp[(i, j)] = RQRp[(j, i)]
                end
            end
        end
    end

    # Create P_t dictionary
    P_dict = Dict{Tuple{Int64, Int64}, Matrix{S}}()
    T_size = size(T, 1)
    P_dict[(1, 1)] = P_0[1:T_size, 1:T_size]
    if n_simulate_states > 1
        for i = 2:n_simulate_states
            for j = i:n_simulate_states
                P_dict[(i, j)] = P_0[1 + T_size * (i-1):T_size * i, 1 + T_size * (j-1):T_size * j]
            end
        end
    end

    return CTBlockKalmanFilter(T_powers, R_blocks, C, Q, Z, D, E, RQRp,
                               s_0, P_dict, n_simulate_states, NaN)
end


"""
```
ct_block_kalman_filter(y, T, R, C, Q, Z, D, E, s_0 = Vector(), P_0 = Matrix();
    outputs = [:loglh, :pred, :filt], Nt0 = 0)
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
- `T`: `Ns` x `Ns` state transition matrix
- `R`: `Ns` x `Ne` matrix in the transition equation mapping shocks to states
- `C`: `Ns` x 1 constant vector in the transition equation
- `Q`: `Ne` x `Ne` matrix of shock covariances
- `Z`: `Ny` x `Ns` matrix in the measurement equation mapping states to
  observables
- `D`: `Ny` x 1 constant vector in the measurement equation
- `E`: `Ny` x `Ny` matrix of measurement error covariances

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

Compute the one-step-ahead states s_{t|t-1} and state covariances P_{t|t-1} and
assign to `k`.
"""
function ct_block_kalman_filter(y::Matrix{S}, T::Matrix{S}, R::Matrix{S}, C::Vector{S},
                          Q::Matrix{S}, Z::Matrix{S}, D::Vector{S}, E::Matrix{S};
                          n_simulate_states::Int = 1, s_0::Vector{S} = Vector{S}(undef, 0),
                          P_0::Matrix{S} = Matrix{S}(undef, 0, 0),
                          outputs::Vector{Symbol} = [:loglh, :pred, :filt],
                          Nt0::Int = 0) where {S<:AbstractFloat}
    # Check s_0, P_0 right dimensions
    if size(s_0, 1) != n_simulate_states * size(T, 1)
        error("Dimension error: s_0 does not stack intermediate simulated states.")
    elseif size(P_0, 1) != n_simulate_states * size(T, 1)
        error("Dimension error: P_0 does correspond to a state vector with stacked intermediate simulated states.")
    end

    # Determine outputs
    return_loglh = :loglh in outputs
    return_pred  = :pred  in outputs
    return_filt  = :filt  in outputs

    # Initialize inputs and outputs
    k = CTBlockKalmanFilter(T, R, C, Q, Z, D, E, s_0, P_0; n_simulate_states)

    # Dimensions
    Ns = size(T, 1) * n_simulate_states # number of states
    Nt = size(y, 2)                     # number of periods of data

    mynan = convert(S, NaN)
    loglh  = return_loglh ? fill(mynan, Nt)         : Vector{S}(0)
    s_pred = return_pred  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_pred = return_pred  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)
    s_filt = return_filt  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_filt = return_filt  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)

    # Populate initial states
    s_0 = k.s_t
    P_0 = k.P_t # note that this is now a dict

    # Loop through periods t
    for t = 1:Nt
        # Forecast
        forecast!(k)
        if return_pred
            s_pred[:,    t] = k.s_t
            P_tmp = zeros(Ns, Ns)
            T_size = size(T, 1)
            for i = 1:n_simulate_states
                for j = 1:n_simulate_states
                    if i < j
                        P_tmp[1 + T_size * (i-1):T_size * i, 1 + T_size * (j-1):T_size * j] = k.P_t[(j,i)]'
                    else
                        P_tmp[1 + T_size * (i-1):T_size * i, 1 + T_size * (j-1):T_size * j] = k.P_t[(i,j)]
                    end
                end
            end
            P_pred[:, :, t] = P_tmp
        end

        # Update and compute log-likelihood
        update!(k, y[:, t]; return_loglh = return_loglh)
        if return_filt
            s_filt[:, t] = k.s_t
            P_tmp = zeros(Ns, Ns)
            T_size = size(T, 1)
            for i = 1:n_simulate_states
                for j = 1:n_simulate_states
                    if i < j # get transpose of upper triangle
                        P_tmp[1 + T_size * (i-1):T_size * i, 1 + T_size * (j-1):T_size * j] = k.P_t[(j,i)]'
                    else     # held in P_t as an element
                        P_tmp[1 + T_size * (i-1):T_size * i, 1 + T_size * (j-1):T_size * j] = k.P_t[(i,j)]
                    end
                end
            end
            P_filt[:, :, t] = P_tmp
        end
        if return_loglh
            loglh[t]        = k.loglh_t
        end

        # # Update s_0 and P_0 if Nt0 > 0
        # if t == Nt0
        #     s_0 = k.s_t
        #     P_0 = k.P_t
        # end
    end

    # Populate final states
    s_T = k.s_t
    P_T = k.P_t

    # Remove presample periods
    loglh, s_pred, P_pred, s_filt, P_filt =
        remove_presample!(Nt0, loglh, s_pred, P_pred, s_filt, P_filt; outputs = outputs)

    return loglh, s_pred, P_pred, s_filt, P_filt, s_0, P_0, s_T, P_T
end

# Computes the one-step-ahead states s_{t|t-1} and state covariances P_{t|t-1}
# and assign to 'k'
function forecast!(k::CTBlockKalmanFilter{Float64})
    T_powers, C, RQRp = k.T_powers, k.C, k.RQRp
    orig_dim = size(T_powers[1], 1)
    s_filt_prev = copy(k.s_t[end - orig_dim:orig_dim]) # don't need all of s_filt
    P_filt = k.P_t

    for i = 1:k.n_simulate_states
        k.s_t[1 + (i - 1) * orig_dim:orig_dim * i] = T_powers[i] * s_filt_prev + C
        for j = i:k.n_simulate_states
            k.P_t[(i,j)] = T_powers[i] * P_filt[(n_simulate_states, n_simulate_states)] * T_powers[j]' + RQRp[(i,j)]
        end
    end

    return nothing
end

function update!(k::CTBlockKalmanFilter{Float64}, y_obs::Vector{Float64}; return_loglh::Bool = true)
    # Keep rows of measurement equation corresponding to non-NaN observables
    nonnan = .!isnan.(y_obs)
    y_obs = y_obs[nonnan]
    Z = k.Z[nonnan, :]
    D = k.D[nonnan]
    E = k.E[nonnan, nonnan]
    Ny = length(y_obs)
    s_pred = k.s_t
    P_pred = k.P_t
    T_size = size(k.T_powers[1], 1)

    y_pred = Z * s_pred  + D
    V_pred = Z * P_pred * Z' + E
    V_pred = (V_pred + V_pred')/2
    V_pred_inv = inv(V_pred)
    dy = y_obs - y_pred
    PZV = P_pred' * Z' * V_pred_inv

    k.s_t = s_pred * PZV * dy
    k.P_t = P_pred - PZV * Z * P_pred

    if return_loglh
        k.loglh_t = -(Ny * log(2*π) + log(det(V_pred)) + dy'*V_pred_inv*dy)/2
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
