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
struct Kalman{S<:AbstractFloat}
    loglh::Vector{S}            # log p(y_t | y_{1:t-1}), t = 1:T
    s_pred::AbstractArray{S}    # s_{t|t-1}, t = 1:T
    P_pred::Array{S, 3}         # P_{t|t-1}, t = 1:T
    s_filt::AbstractArray{S}    # s_{t|t}, t = 1:T
    P_filt::Array{S, 3}         # P_{t|t}, t = 1:T
    s_0::Vector{S}              # s_0
    P_0::AbstractArray{S}       # P_0
    s_T::Vector{S}              # s_{T|T}
    P_T::AbstractArray{S}       # P_{T|T}
    total_loglh::S              # log p(y_{1:t})
end

function Kalman(loglh::Vector{S},
                s_pred::AbstractArray{S}, P_pred::Array{S, 3},
                s_filt::AbstractArray{S}, P_filt::Array{S, 3},
                s_0::Vector{S}, P_0::AbstractArray{S},
                s_T::Vector{S}, P_T::AbstractArray{S}) where S<:AbstractFloat

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

    return DSGE.Kalman(kal[:loglh][inds],        # loglh
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

function Base.cat(m::AbstractDSGEModel, k1::Kalman{S},
                  k2::Kalman{S}; allout::Bool = true) where S<:AbstractFloat

    loglh  = cat(k1[:loglh], k2[:loglh], dims = 1)
    s_pred = cat(k1[:s_pred], k2[:s_pred], dims = 2)
    P_pred = cat(k1[:P_pred], k2[:P_pred], dims = 3)
    s_filt = cat(k1[:s_filt], k2[:s_filt], dims = 2)
    P_filt = cat(k1[:P_filt], k2[:P_filt], dims = 3)
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

Returns a Vector{UnitRange{Int64}} of index ranges for the pre- and post-ZLB
regimes. The optional argument `start_date` indicates the first quarter of
`data`.
"""
function zlb_regime_indices(m::AbstractDSGEModel{S}, data::AbstractArray,
                            start_date::Dates.Date=date_presample_start(m)) where S<:AbstractFloat
    T = size(data, 2)
    if n_mon_anticipated_shocks(m) > 0 && !isempty(data)
        if start_date < date_presample_start(m)
            error("Start date $start_date must be >= date_presample_start(m)")

        elseif 0 < subtract_quarters(date_zlb_start(m), start_date) < T
            n_nozlb_periods = subtract_quarters(date_zlb_start(m), start_date)
            regime_inds::Vector{UnitRange{Int64}} = [1:n_nozlb_periods, (n_nozlb_periods+1):T]
        else
            regime_inds = UnitRange{Int64}[1:T]
        end
    else
        regime_inds = UnitRange{Int64}[1:T]
    end
    return regime_inds
end

"""
```
function zlb_plus_regime_indices(m::AbstractDSGEModel{S}, data::AbstractArray,
                                 start_date::Dates.Date=date_presample_start(m)) where S<:AbstractFloat
```
returns (1) a Vector{UnitRange{Int64}} of index ranges for regime switches with the pre- and post-ZLB
regimes spliced into the regime switches and (2) the index in the vector indicating where the ZLB starts.
The optional argument `start_date` indicates the first quarter of
`data`. Use the `Setting` with key `:regime_dates` to set the start dates of different regimes (excluding
the ZLB regime), and use the `Setting` key `:date_zlb_start` to set the start of the post-ZLB regime.
"""
function zlb_plus_regime_indices(m::AbstractDSGEModel{S}, data::AbstractArray,
                                 start_date::Dates.Date=date_presample_start(m)) where S<:AbstractFloat

    T = size(data, 2)
    if !isempty(data)
        # Calculate the number of periods since start date for each regime
        n_regime_periods = Vector{Int}(undef, length(get_setting(m, :regime_dates)))
        for (k, v) in get_setting(m, :regime_dates)
            n_regime_periods[k] = subtract_quarters(v, start_date)
        end

        if start_date < date_presample_start(m)
            error("Start date $start_date must be >= date_presample_start(m)")
        elseif 0 < subtract_quarters(date_zlb_start(m), start_date) < T
            regime_inds = Vector{UnitRange{Int64}}(undef, length(n_regime_periods) + 1)

            n_nozlb_periods  = subtract_quarters(date_zlb_start(m), start_date) # number of periods since start date for start of ZLB

            # Get index of next regime after ZLB starts.
            # Note that it cannot be 1 b/c the first regime starts at the start date
            i_splice_zlb = findfirst(n_nozlb_periods .< n_regime_periods)

            # Populate vector of regime indices
            if isnothing(i_splice_zlb) # post-ZLB is the last regime
                for reg in 1:(length(n_regime_periods) - 1)
                    regime_inds[reg] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
                end
                regime_inds[end - 1] = (n_regime_periods[end] + 1):n_nozlb_periods
                regime_inds[end]     = (n_nozlb_periods + 1):T

                i_zlb_start = length(regime_inds)
            elseif i_splice_zlb > 2 # at least one full regime before ZLB starts
                regime_inds[1] = 1:n_regime_periods[2]
                for reg in 2:i_splice_zlb - 2 # if i_splice_zlb == 3, then this loop does not run
                    regime_inds[reg] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
                end
                regime_inds[i_splice_zlb - 1] = (n_regime_periods[i_splice_zlb - 2] + 1):n_nozlb_periods
                regime_inds[i_splice_zlb]     = (n_nozlb_periods + 1):(n_regime_periods[i_splice_zlb - 1])

                i_zlb_start = i_splice_zlb

                # Index reg + 1 b/c we have spliced pre- and post- ZLB regime in
                for reg in i_splice_zlb:(length(n_regime_periods) - 1)
                    regime_inds[reg + 1] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
                end
                regime_inds[end] = (n_regime_periods[end] + 1):T
            else # post- ZLB regimes start in the very first regime specified by user
                regime_inds[1] = 1:n_nozlb_periods
                regime_inds[2] = (n_nozlb_periods + 1):n_regime_periods[1]

                i_zlb_start = 2

                # Index reg + 1 b/c spliced pre- and post-ZLB regime in
                for reg in 2:length(n_regime_periods)
                    regime_inds[reg + 1] = (n_regime_periods[reg - 1] + 1):n_regime_periods[reg]
                end
            end
        else
            # This is the case that date_zlb_start <= start_date so the first regime is the post-ZLB regime (no pre-ZLB)
            # Since there is no ZLB regime to splice in, we return the regime_inds corresponding to the user-specified regimes
            regime_inds = Vector{UnitRange{Int64}}(undef, length(n_regime_periods))
            for reg in 1:(length(n_regime_periods) - 1)
                regime_inds[reg] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
            end
            regime_inds[end] = (n_regime_periods[end] + 1):T

            i_zlb_start = 1 # start immediately in post-ZLB
        end
    else # Empty, so we ignore regime switching
        regime_inds = UnitRange{Int64}[1:T]
    end
    return regime_inds, i_zlb_start
end

"""
```
function regime_indices(m::AbstractDSGEModel{S}, start_date::Dates.Date = date_presample_start(m),
                        end_date::Dates.Date = prev_quarter(date_forecast_start(m))) where S<:AbstractFloat
```
returns a `Vector{UnitRange{Int64}}` of index ranges for regime switches
The optional argument `start_date` indicates the first quarter of
`data`. The optimal argument `end_date` indicates the last quarter of `data`.
Use the `Setting` with key `:regime_dates` to set the start dates of different regimes.

This function differs from `zlb_plus_regime_indices` because it does not splice the
pre- and post-ZLB regimes into the output of index ranges for regime switches. When it is not necessary
to explicitly zero out anticipated monetary policy shocks to establish the pre-ZLB regime, there is no need
to distinguish between the pre-ZLB and post-ZLB regimes. For example, with shock decompositions,
the user needs to specify the history of filtered shocks, which should have zeroes already for anticipated shocks
in the pre-ZLB regime. Thus, there is no need to explicitly account for the pre-ZLB regime.
"""
function regime_indices(m::AbstractDSGEModel{S}, start_date::Dates.Date = date_presample_start(m),
                        end_date::Dates.Date = prev_quarter(date_forecast_start(m))) where S<:AbstractFloat
    T = subtract_quarters(end_date, start_date) + 1

    n_regime_periods = Vector{Int}(undef, length(get_setting(m, :regime_dates)))
    for (k, v) in get_setting(m, :regime_dates)
        n_regime_periods[k] = subtract_quarters(v, start_date)
    end
    regime_inds = Vector{UnitRange{Int64}}(undef, length(n_regime_periods))

    for reg in 1:(length(n_regime_periods) - 1)
        regime_inds[reg] = (n_regime_periods[reg] + 1):n_regime_periods[reg + 1]
    end
    regime_inds[end] = (n_regime_periods[end] + 1):T

    return regime_inds
end

"""
```
zlb_regime_matrices(m, system, start_date = date_presample_start(m))
```
Returns `TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs`, an 8-tuple of
`Vector{AbstractArray{S}}`s and `Vector{Vector{S}}`s of system matrices for the pre-
and post-ZLB regimes. Of these, only `QQ` changes from pre- to post-ZLB: the
entries corresponding to anticipated shock variances are zeroed out pre-ZLB.
"""
function zlb_regime_matrices(m::AbstractDSGEModel{S}, system::System{S},
                             start_date::Dates.Date=date_presample_start(m)) where S<:AbstractFloat
    if n_mon_anticipated_shocks(m) > 0
        if start_date < date_presample_start(m)
            error("Start date $start_date must be >= date_presample_start(m)")

        # TODO: This technically doesn't handle the case where the end_date of the sample
        # is before the start of the ZLB
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

function zlb_plus_regime_matrices(m::AbstractDSGEModel{S}, system::RegimeSwitchingSystem{S},
                                  start_date::Dates.Date = date_presample_start(m);
                                  n_regimes::Int = 0, ind_zlb_start::Int = 0) where S<:AbstractFloat

    if n_regimes <= 0 || ind_zlb_start <= 0
        throw(DomainError())
    end

    if start_date < date_presample_start(m)
        error("Start date $start_date must be >= date_presample_start(m)")

        # TODO: This technically doesn't handle the case where the end_date of the sample
        # is before the start of the ZLB
    elseif date_presample_start(m) <= start_date <= date_zlb_start(m)
        TTTs = Vector{Matrix{S}}(undef, n_regimes)
        RRRs = Vector{Matrix{S}}(undef, n_regimes)
        CCCs = Vector{Vector{S}}(undef, n_regimes)
        ZZs  = Vector{Matrix{S}}(undef, n_regimes)
        DDs  = Vector{Vector{S}}(undef, n_regimes)
        QQs  = Vector{Matrix{S}}(undef, n_regimes)
        EEs  = Vector{Matrix{S}}(undef, n_regimes)

        shock_inds = inds_shocks_no_ant(m)
        if ind_zlb_start == 2
            if n_mon_anticipated_shocks(m) == 0
                QQ_preZLB = copy(system[1, :QQ]) # anticipated shocks are zero already
            else
                QQ_preZLB                         = zeros(size(system[1, :QQ]))
                QQ_preZLB[shock_inds, shock_inds] = system[1, :QQ][shock_inds, shock_inds]
            end
            QQ_postZLB = copy(system[1, :QQ])
            TTTs[1:2] = [system[1, :TTT], system[1, :TTT]]::Vector{Matrix{S}}
            RRRs[1:2] = [system[1, :RRR], system[1, :RRR]]::Vector{Matrix{S}}
            CCCs[1:2] = [system[1, :CCC], system[1, :CCC]]::Vector{Vector{S}}
            ZZs[1:2]  = [system[1, :ZZ],  system[1, :ZZ]]::Vector{Matrix{S}}
            DDs[1:2]  = [system[1, :DD],  system[1, :DD]]::Vector{Vector{S}}
            EEs[1:2]  = [system[1, :EE],  system[1, :EE]]::Vector{Matrix{S}}
            QQs[1]    = QQ_preZLB
            QQs[2]    = QQ_postZLB

            input = 2:n_regimes
            TTTs[3:end] = map(i -> system[i, :TTT], input)::Vector{Matrix{S}}
            RRRs[3:end] = map(i -> system[i, :RRR], input)::Vector{Matrix{S}}
            CCCs[3:end] = map(i -> system[i, :CCC], input)::Vector{Vector{S}}
            ZZs[3:end]  = map(i -> system[i, :ZZ],  input)::Vector{Matrix{S}}
            QQs[3:end]  = map(i -> system[i, :QQ],  input)::Vector{Matrix{S}}
            DDs[3:end]  = map(i -> system[i, :DD],  input)::Vector{Vector{S}}
            EEs[3:end]  = map(i -> system[i, :EE],  input)::Vector{Matrix{S}}
        elseif ind_zlb_start > 1
            # Populate matrices for regimes before the ZLB split
            pre_zlb_input = 1:(ind_zlb_start - 2) # minus 2 b/c ind_zlb_start - 1 will be the pre-ZLB regime
            splice_reg    = ind_zlb_start - 1

            TTTs[pre_zlb_input] = map(i -> system[i, :TTT], pre_zlb_input)::Vector{Matrix{S}}
            RRRs[pre_zlb_input] = map(i -> system[i, :RRR], pre_zlb_input)::Vector{Matrix{S}}
            CCCs[pre_zlb_input] = map(i -> system[i, :CCC], pre_zlb_input)::Vector{Vector{S}}
            ZZs[pre_zlb_input]  = map(i -> system[i, :ZZ],  pre_zlb_input)::Vector{Matrix{S}}
            DDs[pre_zlb_input]  = map(i -> system[i, :DD],  pre_zlb_input)::Vector{Vector{S}}
            EEs[pre_zlb_input]  = map(i -> system[i, :EE],  pre_zlb_input)::Vector{Matrix{S}}

            if n_mon_anticipated_shocks(m) == 0
                QQs[pre_zlb_input] = map(i -> system[i, :QQ], pre_zlb_input)::Vector{Matrix{S}}
            else
                for i in pre_zlb_input
                    QQ_preZLB                         = zeros(size(system[i, :QQ]))
                    QQ_preZLB[shock_inds, shock_inds] = system[i, :QQ][shock_inds, shock_inds]
                    QQs[i]                            = QQ_preZLB
                end
            end

            # Regime in which the ZLB occurs
            if n_mon_anticipated_shocks(m) == 0
                QQ_preZLB = copy(system[splice_reg, :QQ]) # anticipated shocks are zero already
            else
                QQ_preZLB                         = zeros(size(system[splice_reg, :QQ]))
                QQ_preZLB[shock_inds, shock_inds] = system[splice_reg, :QQ][shock_inds, shock_inds]
            end
            QQ_postZLB = copy(system[splice_reg, :QQ])
            TTTs[splice_reg:ind_zlb_start] = [system[splice_reg, :TTT], system[splice_reg, :TTT]]::Vector{Matrix{S}}
            RRRs[splice_reg:ind_zlb_start] = [system[splice_reg, :RRR], system[splice_reg, :RRR]]::Vector{Matrix{S}}
            CCCs[splice_reg:ind_zlb_start] = [system[splice_reg, :CCC], system[splice_reg, :CCC]]::Vector{Vector{S}}
            ZZs[splice_reg:ind_zlb_start]  = [system[splice_reg, :ZZ],  system[splice_reg, :ZZ]]::Vector{Matrix{S}}
            DDs[splice_reg:ind_zlb_start]  = [system[splice_reg, :DD],  system[splice_reg, :DD]]::Vector{Vector{S}}
            EEs[splice_reg:ind_zlb_start]  = [system[splice_reg, :EE],  system[splice_reg, :EE]]::Vector{Matrix{S}}
            QQs[splice_reg]     = QQ_preZLB
            QQs[splice_reg + 1] = QQ_postZLB

            # Handle regimes after the ZLB split
            post_zlb_input = (ind_zlb_start + 1):n_regimes

            TTTs[post_zlb_input] = map(i -> system[i, :TTT], post_zlb_input)::Vector{Matrix{S}}
            RRRs[post_zlb_input] = map(i -> system[i, :RRR], post_zlb_input)::Vector{Matrix{S}}
            CCCs[post_zlb_input] = map(i -> system[i, :CCC], post_zlb_input)::Vector{Vector{S}}
            ZZs[post_zlb_input]  = map(i -> system[i, :ZZ],  post_zlb_input)::Vector{Matrix{S}}
            QQs[post_zlb_input]  = map(i -> system[i, :QQ],  post_zlb_input)::Vector{Matrix{S}}
            DDs[post_zlb_input]  = map(i -> system[i, :DD],  post_zlb_input)::Vector{Vector{S}}
            EEs[post_zlb_input]  = map(i -> system[i, :EE],  post_zlb_input)::Vector{Matrix{S}}
        end
    elseif date_zlb_start(m) < start_date
        # We always use post-ZLB matrices in this case,
        # hence we just return the matrices of the `RegimeSwitchingSystem`
        input = 1:n_regimes
        TTTs = map(i -> system[i, :TTT], input)::Vector{Matrix{S}}
        RRRs = map(i -> system[i, :RRR], input)::Vector{Matrix{S}}
        CCCs = map(i -> system[i, :CCC], input)::Vector{Vector{S}}
        ZZs  = map(i -> system[i, :ZZ],  input)::Vector{Matrix{S}}
        QQs  = map(i -> system[i, :QQ],  input)::Vector{Matrix{S}}
        DDs  = map(i -> system[i, :DD],  input)::Vector{Vector{S}}
        EEs  = map(i -> system[i, :EE],  input)::Vector{Matrix{S}}
    end


    return TTTs, RRRs, CCCs, QQs, ZZs, DDs, EEs
end
