"""
```
Kalman{S<:AbstractFloat}
```

### Fields:

- `L`: value of the average log likelihood function of the SSM under assumption that
  observation noise ϵ(t) is normally distributed
- `zend`: state vector in the last period for which data is provided
- `Pend`: variance-covariance matrix for `zend`

#### Fields filled in when `allout` in a call to `kalman_filter`:

- `pred`: a `Nz` x `T` matrix containing one-step predicted state vectors
- `vpred`: a `Nz` x `Nz` x `T` matrix containing mean square errors of predicted
  state vectors
- `filt`: an `Nz` x `T` matrix containing filtered state vectors
- `vfilt`: an `Nz` x `Nz` x `T` matrix containing mean square errors of filtered
  state vectors
- `z0`: starting-period state vector. If there are presample periods in the
  data, then `z0` is the state vector at the end of the presample/beginning of
  the main sample
- `P0`: variance-covariance matrix for `z0`
"""
immutable Kalman{S<:AbstractFloat}
    L::S                  # likelihood
    zend::Vector{S}       # last-period state vector
    Pend::Matrix{S}       # last-period variance-covariance matrix for the states
    pred::Matrix{S}       # predicted value of states in period T+1
    vpred::Array{S, 3}    # predicted variance-covariance matrix for states in period T+1
    yprederror::Matrix{S}
    ystdprederror::Matrix{S}
    rmse::Matrix{S}
    rmsd::Matrix{S}
    filt::Matrix{S}        # filtered states
    vfilt::Array{S, 3}     # mean squared errors of filtered state vectors
    z0::Vector{S}          # starting-period state vector
    vz0::Matrix{S}         # starting-period variance-covariance matrix for the states
end

function Kalman{S<:AbstractFloat}(L::S,
                                  zend::Vector{S}          = Vector{S}(),
                                  Pend::Matrix{S}          = Matrix{S}(),
                                  pred::Matrix{S}          = Matrix{S}(),
                                  vpred::Array{S, 3}       = Array{S}(0, 0, 0),
                                  filt::Matrix{S}          = Matrix{S}(),
                                  vfilt::Array{S, 3}       = Array{S}(0, 0, 0),
                                  yprederror::Matrix{S}    = Matrix{S}(),
                                  ystdprederror::Matrix{S} = Matrix{S}(),
                                  rmse::Matrix{S}          = Matrix{S}(),
                                  rmsd::Matrix{S}          = Matrix{S}(),
                                  z0::Vector{S}            = Vector{S}(),
                                  P0::Matrix{S}            = Matrix{S}())

    return Kalman{S}(L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt, z0, P0)
end

function Base.getindex(K::Kalman, d::Symbol)
    if d in (:L, :zend, :Pend, :pred, :vpred, :yprederror, :ystdprederror, :rmse, :rmsd,
             :filt, :vfilt, :z0, :vz0)
        return getfield(K, d)
    else
        throw(KeyError(d))
    end
end

function Base.cat{S<:AbstractFloat}(m::AbstractModel, k1::Kalman{S},
    k2::Kalman{S}; allout::Bool = true)

    L = k1[:L] + k2[:L]
    zend = k2[:zend]
    Pend = k2[:Pend]

    if allout
        pred  = hcat(k1[:pred], k2[:pred])
        vpred = cat(3, k1[:vpred], k2[:vpred])
        yprederror    = hcat(k1[:yprederror], k2[:yprederror])
        ystdprederror = hcat(k1[:ystdprederror], k2[:yprederror])
        rmse  = sqrt(mean((yprederror.^2)', 1))
        rmsd  = sqrt(mean((ystdprederror.^2)', 1))
        filt  = hcat(k1[:filt], k2[:filt])
        vfilt = cat(3, k1[:vfilt], k2[:vfilt])
        z0    = k1[:z0]
        P0    = k1[:vz0]

        return Kalman(L, zend, Pend, pred, vpred, yprederror, ystdprederror,
            rmse, rmsd, filt, vfilt, z0, P0)
    else
        return Kalman(L, zend, Pend)
    end
end

"""
```
zlb_regime_indices(m, data)
```

Returns a Vector{Range{Int64}} of index ranges for the pre- and post-ZLB regimes.
"""
function zlb_regime_indices{S<:AbstractFloat}(m::AbstractModel{S}, data::Matrix{S})
    T = size(data, 2)

    if n_anticipated_shocks(m) > 0
        regime_inds = Vector{Range{Int64}}(2)
        regime_inds[1] = 1:index_zlb_start(m)-1
        regime_inds[2] = index_zlb_start(m):T # allows for conditional data
    else
        regime_inds = Range{Int64}[1:T]
    end

    return regime_inds
end

"""
```
zlb_regime_matrices(m, system)
```

Returns `TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs`, an 8-tuple of
`Vector{Matrix{S}}`s and `Vector{Vector{S}}`s of system matrices for the pre-
and post-ZLB regimes. Of these, only `QQ` changes from pre- to post-ZLB: the
entries corresponding to anticipated shock variances are zeroed out pre-ZLB.
"""
function zlb_regime_matrices{S<:AbstractFloat}(m::AbstractModel{S}, system::System{S})
    if !all(x -> x == 0, system[:MM])
        error("Kalman filter and smoothers not implemented for nonzero MM")
    end

    if n_anticipated_shocks(m) > 0
        n_regimes = 2

        shock_inds = inds_shocks_no_ant(m)
        QQ_ZLB = system[:QQ]
        QQ_preZLB = zeros(size(QQ_ZLB))
        QQ_preZLB[shock_inds, shock_inds] = QQ_ZLB[shock_inds, shock_inds]
        QQs = Matrix{S}[QQ_preZLB, QQ_ZLB]
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