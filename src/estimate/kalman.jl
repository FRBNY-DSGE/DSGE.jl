"""
```
Kalman{S<:AbstractFloat}
```
### Fields:

- `loglh`: vector of conditional log-likelihoods log p(y_t | y_{1:t-1}), t = 1:T
- `s_T`: state vector in the last period for which data is provided
- `P_T`: variance-covariance matrix for `s_T`
- `s_pred`: `Ns` x `Nt` matrix of s_{t|t-1}, t = 1:T
- `P_pred`: `Ns` x `Ns` x `Nt` array of P_{t|t-1}, t = 1:T
- `s_filt`: `Ns` x `Nt` matrix of s_{t|t}, t = 1:T
- `P_filt`: `Ns` x `Ns` x `Nt` array of P_{t|t}, t = 1:T
- `s_0`: starting-period state vector. If there are presample periods in the
  data, then `s_0` is the state vector at the end of the presample/beginning of
  the main sample
- `P_0`: variance-covariance matrix for `s_0`
- `total_loglh`: log p(y_{1:t})
"""
immutable Kalman{S<:AbstractFloat}
    loglh::Vector{S}     # log p(y_t | y_{1:t-1}), t = 1:T
    s_pred::Matrix{S}    # s_{t|t-1}, t = 1:T
    P_pred::Array{S, 3}  # P_{t|t-1}, t = 1:T
    s_filt::Matrix{S}    # s_{t|t}, t = 1:T
    P_filt::Array{S, 3}  # P_{t|t}, t = 1:T
    s_0::Vector{S}       # s_0
    P_0::Matrix{S}       # P_0
    s_T::Vector{S}       # s_{T|T}
    P_T::Matrix{S}       # P_{T|T}
    total_loglh::S       # log p(y_{1:t})
end

function Kalman(loglh::Vector{S},
                s_pred::Matrix{S}, P_pred::Array{S, 3},
                s_filt::Matrix{S}, P_filt::Array{S, 3},
                s_0::Vector{S}, P_0::Matrix{S},
                s_T::Vector{S}, P_T::Matrix{S}) where {S<:AbstractFloat}

    return Kalman{S}(loglh, s_pred, P_pred, s_filt, P_filt, s_0, P_0, s_T, P_T, sum(loglh))
end

function Base.getindex(K::Kalman, d::Symbol)
    if d in (:loglh, :s_pred, :P_pred, :s_filt, :P_filt, :s_0, :P_0, :s_T, :P_T, :total_loglh)
        return getfield(K, d)
    else
        throw(KeyError(d))
    end
end

function Base.getindex(kal::Kalman, inds::Union{Int, UnitRange{Int}})
    t0 = first(inds)
    t1 = last(inds)

    return Kalman(kal[:loglh][inds],        # loglh
                  kal[:s_pred][:,    inds], # s_pred
                  kal[:P_pred][:, :, inds], # P_pred
                  kal[:s_filt][:,    inds], # filt
                  kal[:P_filt][:, :, inds], # P_filt
                  kal[:s_filt][:,    t0],   # s_0
                  kal[:P_filt][:, :, t0],   # P_0
                  kal[:s_filt][:,    t1],   # s_T
                  kal[:P_filt][:, :, t1],   # P_T
                  sum(kal[:loglh][inds]))   # total_loglh
end

function Base.cat{S<:AbstractFloat}(m::AbstractModel, k1::Kalman{S},
    k2::Kalman{S}; allout::Bool = true)

    loglh  = cat(1, k1[:loglh], k2[:loglh])
    s_pred = cat(2, k1[:s_pred], k2[:s_pred])
    P_pred = cat(3, k1[:P_pred], k2[:P_pred])
    s_filt = cat(2, k1[:s_filt], k2[:s_filt])
    P_filt = cat(3, k1[:P_filt], k2[:P_filt])
    s_0    = k1[:s_0]
    P_0    = k1[:P_0]
    s_T    = k2[:s_T]
    P_T    = k2[:P_T]
    total_loglh = k1[:total_loglh] + k2[:total_loglh]

    return Kalman(loglh, s_pred, P_pred, s_filt, P_filt, s_0, P_0, s_T, P_T, total_loglh)
end

"""
```
zlb_regime_indices(m, data, start_date = date_presample_start(m))
```

Returns a Vector{Range{Int64}} of index ranges for the pre- and post-ZLB
regimes. The optional argument `start_date` indicates the first quarter of
`data`.
"""
function zlb_regime_indices{S<:AbstractFloat}(m::AbstractModel{S}, data::Matrix{S},
                                              start_date::Date = date_presample_start(m))
    T = size(data, 2)

    if n_anticipated_shocks(m) > 0 && !isempty(data)
        if start_date < date_presample_start(m)
            error("Start date $start_date must be >= date_presample_start(m)")

        elseif date_presample_start(m) <= start_date <= date_zlb_start(m)
            n_nozlb_periods = subtract_quarters(date_zlb_start(m), start_date)
            regime_inds = Vector{Range{Int64}}(2)
            regime_inds[1] = 1:n_nozlb_periods
            regime_inds[2] = (n_nozlb_periods+1):T

        else # date_zlb_start(m) < start_date
            regime_inds = Range{Int64}[1:T]
        end
    else
        regime_inds = Range{Int64}[1:T]
    end

    return regime_inds
end

"""
```
zlb_regime_matrices(m, system, start_date = date_presample_start(m))
```
Returns `TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs`, an 8-tuple of
`Vector{Matrix{S}}`s and `Vector{Vector{S}}`s of system matrices for the pre-
and post-ZLB regimes. Of these, only `QQ` changes from pre- to post-ZLB: the
entries corresponding to anticipated shock variances are zeroed out pre-ZLB.
"""
function zlb_regime_matrices{S<:AbstractFloat}(m::AbstractModel{S}, system::System{S},
                                               start_date::Date = date_presample_start(m))
    if n_anticipated_shocks(m) > 0
        if start_date < date_presample_start(m)
            error("Start date $start_date must be >= date_presample_start(m)")

        elseif date_presample_start(m) <= start_date <= date_zlb_start(m)
            n_regimes = 2

            shock_inds = inds_shocks_no_ant(m)
            QQ_ZLB = system[:QQ]
            QQ_preZLB = zeros(size(QQ_ZLB))
            QQ_preZLB[shock_inds, shock_inds] = QQ_ZLB[shock_inds, shock_inds]
            QQs = Matrix{S}[QQ_preZLB, QQ_ZLB]

        elseif date_zlb_start(m) < start_date
            n_regimes = 1
            QQs = Matrix{S}[system[:QQ]]
        end
    else
        n_regimes = 1
        QQs = Matrix{S}[system[:QQ]]
    end

    TTTs = fill(system[:TTT], n_regimes)
    RRRs = fill(system[:RRR], n_regimes)
    CCCs = fill(system[:CCC], n_regimes)
    ZZs  = fill(system[:ZZ], n_regimes)
    DDs  = fill(system[:DD], n_regimes)
    EEs  = fill(system[:EE], n_regimes)

    return TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs
end
