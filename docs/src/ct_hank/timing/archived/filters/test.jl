

mutable struct KalmanFilter{S<:AbstractFloat}
    T::Matrix{S}
    R::Matrix{S}
    C::Vector{S}
    Q::Matrix{S}
    Z::Matrix{S}
    D::Vector{S}
    E::Matrix{S}
    s_t::Vector{S} # s_{t|t-1} or s_{t|t}
    P_t::Matrix{S} # P_{t|t-1} or P_{t|t}
    loglh_t::S     # P(y_t | y_{1:t})
end


function compute_values(y::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Vector{S},
    Q::Matrix{S}, Z::Matrix{S}, D::Vector{S}, E::Matrix{S},
    s_0::Vector{S} = Vector{S}(0), P_0::Matrix{S} = Matrix{S}(0, 0);
    outputs::Vector{Symbol} = [:loglh, :pred, :filt],
    Nt0::Int = 0) where {S<:AbstractFloat}
    #Number of time periods for which we have data
    Nt = size(y,2)
    println(size(Nt))
    true_ll, ~, ~, ~, P_filt, ~, ~, ~, ~ = kalman_filter(Nt,y,T,R,C,Q,Z,D,E,s_0,P_0)
    truelik = sum(true_ll)

    norm_P_T = zeros(Nt)
    ch_ll = zeros(Nt)
#    loglh = zeros(Nt)

    #temp var tracks Pt-1|t-2
    #nt is the t that we start using barT
    temp = Array{Float64}(size(P_0)...)
    P_T = 0
    for nt in 1:Nt
        loglh, ~, ~, ~, ~, ~, ~, ~, P_T = kalman_filter(nt,y,T,R,C,Q,Z,D,E,s_0,P_0)
        if nt == 1
            temp = copy(P_T)
        else
            lik = sum(loglh)
            #compute norm
            norm_P_T[nt] = norm(P_T .- temp,2)
            #update temp
            temp = copy(P_T)
            #compute change in likelihood
            ch_ll[nt] = abs(truelik - lik)
        end
    end
    return norm_P_T, ch_ll, truelik
end

#updated inputs to track nt
function kalman_filter(nt::Int64, y::Matrix{S},
    T::Matrix{S}, R::Matrix{S}, C::Vector{S},
    Q::Matrix{S}, Z::Matrix{S}, D::Vector{S}, E::Matrix{S},
    s_0::Vector{S} = Vector{S}(0), P_0::Matrix{S} = Matrix{S}(0, 0);
    outputs::Vector{Symbol} = [:loglh, :pred, :filt],
    Nt0::Int = 0) where {S<:AbstractFloat}

    # Determine outputs
    return_loglh = :loglh in outputs
    return_pred  = :pred  in outputs
    return_filt  = :filt  in outputs

    # Dimensions
    Ns = size(T, 1) # number of states
    Nt = size(y, 2) # number of periods of data

    # Initialize inputs and outputs
    k = KalmanFilter(T, R, C, Q, Z, D, E, s_0, P_0)

    mynan = convert(S, NaN)
    loglh  = return_loglh ? fill(mynan, Nt)         : Vector{S}(0)
    s_pred = return_pred  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_pred = return_pred  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)
    s_filt = return_filt  ? fill(mynan, Ns, Nt)     : Matrix{S}(0, 0)
    P_filt = return_filt  ? fill(mynan, Ns, Ns, Nt) : Array{S, 3}(0, 0, 0)

    # Populate initial states
    s_0 = k.s_t
    P_0 = k.P_t

    # Loop through periods t
    for t = 1:Nt
        # Forecast
        forecast!(k,t,nt)
        if return_pred
            s_pred[:,    t] = k.s_t
            P_pred[:, :, t] = k.P_t
        end

        # Update and compute log-likelihood

# update inputs to track t, nt
        update!(t, nt, k, y[:, t]; return_loglh = return_loglh)
        if return_filt
            s_filt[:,    t] = k.s_t
            P_filt[:, :, t] = k.P_t
        end
        if return_loglh
            loglh[t]        = k.loglh_t
        end

        # Update s_0 and P_0 if Nt0 > 0
        if t == Nt0
            s_0 = k.s_t
            P_0 = k.P_t
        end
    end

    # Populate final states
    s_T = k.s_t
    P_T = k.P_t

    # Remove presample periods
    loglh, s_pred, P_pred, s_filt, P_filt =
        remove_presample!(Nt0, loglh, s_pred, P_pred, s_filt, P_filt; outputs = outputs)

    return loglh, s_pred, P_pred, s_filt, P_filt, s_0, P_0, s_T, P_T
end

function forecast!(k::KalmanFilter{S}, t::Int64, nt::Int64) where {S<:AbstractFloat}
    T, R, C, Q = k.T, k.R, k.C, k.Q
    s_filt, P_filt = k.s_t, k.P_t

    k.s_t = T*s_filt + C         # s_{t|t-1} = T*s_{t-1|t-1} + C
    if t < nt
        k.P_t = T*P_filt*T' + R*Q*R' # P_{t|t-1} = Var s_{t|t-1} = T*P_{t-1|t-1}*T' + R*Q*R'
    end
    return nothing
end

#update inputs
function update!(t::Int64, nt::Int64, k::KalmanFilter{S}, y_obs::Vector{S};
                 return_loglh::Bool = true) where {S<:AbstractFloat}
    # Keep rows of measurement equation corresponding to non-NaN observables
    nonnan = .!isnan.(y_obs)
    y_obs = y_obs[nonnan]
    Z = k.Z[nonnan, :]
    D = k.D[nonnan]
    E = k.E[nonnan, nonnan]
    Ny = length(y_obs)

    s_pred = k.s_t
    P_pred = k.P_t

    y_pred = Z*s_pred + D         # y_{t|t-1} = Z*s_{t|t-1} + D
    V_pred = Z*P_pred*Z' + E      # V_{t|t-1} = Var y_{t|t-1} = Z*P_{t|t-1}*Z' + E
    V_pred = (V_pred + V_pred')/2
    V_pred_inv = inv(V_pred)
    dy = y_obs - y_pred           # dy = y_t - y_{t|t-1} = prediction error
    PZV = P_pred'*Z'*V_pred_inv

    k.s_t = s_pred + PZV*dy       # s_{t|t} = s_{t|t-1} + P_{t|t-1}'*Z'/V_{t|t-1}*dy

    #if t is less than nt, keep updating P_pred as usual
    if t < nt
        k.P_t = P_pred - PZV*Z*P_pred # P_{t|t} = P_{t|t-1} - P_{t|t-1}'*Z'/V_{t|t-1}*Z*P_{t|t-1}
    #we are using SS value of P_t
    else
        k.P_t = P_pred
    end

    if return_loglh
        k.loglh_t = -(Ny*log(2π) + log(det(V_pred)) + dy'*V_pred_inv*dy)/2 # p(y_t | y_{1:t-1})
    end
    return nothing
end

function remove_presample!(Nt0::Int, loglh::Vector{S},
                           s_pred::Matrix{S}, P_pred::Array{S, 3},
                           s_filt::Matrix{S}, P_filt::Array{S, 3};
                           outputs::Vector{Symbol} = [:loglh, :pred, :filt]) where {S<:AbstractFloat}
    if Nt0 > 0
        Nt = length(loglh)
        insample = (Nt0+1):Nt

        if :loglh in outputs
            loglh  = loglh[insample]
        end
        if :pred in outputs
            s_pred = s_pred[:,    insample]
            P_pred = P_pred[:, :, insample]
        end
        if :filt in outputs
            s_filt = s_filt[:,    insample]
            P_filt = P_filt[:, :, insample]
        end
    end
    return loglh, s_pred, P_pred, s_filt, P_filt
end