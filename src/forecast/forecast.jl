"""
```
forecast(m, system, z0; enforce_zlb = false, shocks = Matrix{S}(undef, 0,0))

forecast(system, z0, shocks; enforce_zlb = false)
```

### Inputs

- `system::System{S}`: state-space system matrices
- `z0::Vector{S}`: state vector in the final historical period
- `shocks::Matrix{S}`: `nshocks` x `nperiods` matrix of shocks to use when
  forecasting. Note that in the first method, `nperiods` doesn't necessarily
  have to equal `forecast_horizons(m)`; it will be truncated or padded with
  zeros appropriately

**Method 1 only:**

- `m::AbstractDSGEModel`

where `S<:AbstractFloat`.

### Keyword Arguments

- `cond_type::Symbol`: one of `:none`, `:semi`, or `:full`, used to determine
  how many periods to forecast ahead. If `cond_type in [:semi, :full]`, the
  forecast horizon is reduced by the number of periods of conditional
  data. Defaults to `:none`.
- `enforce_zlb::Bool`: whether to enforce the zero lower bound. Defaults to
  `false`.
- `shocks::Matrix{S}`: matrix of size `nshocks` x `shock_horizon` of shock
  innovations under which to forecast. If `shock_horizon > horizon`, the extra
  periods of shocks will be ignored; if `shock_horizon < horizon`, zeros will be
  filled in for the shocks hitting the remaining forecasted periods.
- `draw_shocks::Bool`: if `isempty(shocks)`, indicates whether to draw shocks
  according to:

  1. If `forecast_tdist_shocks(m)`, draw `horizons` many shocks from a
     `Distributions.TDist(forecast_tdist_df_val(m))`
  2. Otherwise, draw `horizons` many shocks from a
     `DegenerateMvNormal(zeros(nshocks), sqrt(system[:QQ]))`

  or to set `shocks` to a `nshocks` x `horizon` matrix of zeros. Defaults to
  `false`. If `shocks` is provided as a keyword argument, this flag has no
  effect.

### Outputs

- `states::Matrix{S}`: matrix of size `nstates` x `horizon` of forecasted states
- `obs::Matrix{S}`: matrix of size `nobs` x `horizon` of forecasted observables
- `pseudo::Matrix{S}`: matrix of size `npseudo` x `horizon` of forecasted
  pseudo-observables
- `shocks::Matrix{S}`: matrix of size `nshocks` x `horizon` of shock innovations
"""
function forecast(m::AbstractDSGEModel, system::Union{RegimeSwitchingSystem{S}, System{S}},
    z0::Vector{S}; cond_type::Symbol = :none, enforce_zlb::Bool = false,
    shocks::AbstractMatrix{S} = Matrix{S}(undef, 0, 0), draw_shocks::Bool = false) where {S<:AbstractFloat}

    # Numbers of things
    nshocks = n_shocks_exogenous(m)
    horizon = forecast_horizons(m; cond_type = cond_type)

    if isempty(shocks)
        # Populate shocks matrix
        if draw_shocks
            # If Rgeime Switching System, need to get proper QQ in each forecast period.
            if isa(system, RegimeSwitchingSystem)
                shocks = zeros(nshocks, horizon)
                μ = zeros(nshocks)
                # Forecast regime indices (starting with 1 in 1st fcast period, going up to horizon)
                regime_inds = get_fcast_regime_inds(m, horizon, cond_type)
                # The regime in which forecast starts (indexing from presample_start as reg=1)
                sys_ind = cond_type == :none ? get_setting(m, :reg_forecast_start) :
                    max(get_setting(m, :reg_forecast_start), get_setting(m, :reg_post_conditional_end))
                for ts in regime_inds
                    σ = sqrt.(system[sys_ind, :QQ])
                    dist = if forecast_tdist_shocks(m)
                        # Use t-distributed shocks
                        ν = forecast_tdist_df_val(m)
                        DegenerateDiagMvTDist(μ, σ, ν)
                    else
                        # Use normally distributed shocks
                        DegenerateMvNormal(μ, σ)
                    end
                    shocks[:, ts] = rand(dist, length(ts))
                    sys_ind += 1
                end
            else
                μ = zeros(S, nshocks)
                σ = sqrt.(system[:QQ])
                dist = if forecast_tdist_shocks(m)
                    # Use t-distributed shocks
                    ν = forecast_tdist_df_val(m)
                    DegenerateDiagMvTDist(μ, σ, ν)
                else
                    # Use normally distributed shocks
                    DegenerateMvNormal(μ, σ)
                end
                shocks = rand(dist, horizon)
            end

            # Forecast without anticipated shocks
            if n_mon_anticipated_shocks(m) > 0
                ind_ant1 = m.exogenous_shocks[:rm_shl1]
                ind_antn = m.exogenous_shocks[Symbol("rm_shl$(n_mon_anticipated_shocks(m))")]
                ant_shock_inds = ind_ant1:ind_antn
                shocks[ant_shock_inds, :] .= 0
            end
        else
            shocks = zeros(S, nshocks, horizon)
        end
    else
        # Adjust size of shocks matrix, padding with zeros or cutting off
        # periods of shocks if necessary
        shock_horizon = size(shocks, 2)
        if shock_horizon <= horizon
            shocks0 = zeros(nshocks, horizon - shock_horizon)
            shocks = hcat(shocks, shocks0)
        else
            shocks = shocks[:, 1:horizon]
        end
    end

    if haskey(get_settings(m), :n_periods_no_shocks)
        if get_setting(m, :n_periods_no_shocks) > 0
            shocks[:, 1:get_setting(m, :n_periods_no_shocks)] = zeros(size(shocks, 1), get_setting(m, :n_periods_no_shocks))
        end
    end

    # Populate shocks matrix under alternative policy, if
    # user has specified a function to do so
    alt_policy = alternative_policy(m)

    if (alt_policy.solve != identity &&
        alt_policy.forecast_init != identity)
        shocks, z0 = alt_policy.forecast_init(m, shocks, z0, cond_type = cond_type)
    end

    if haskey(m.settings, :pgap_type) && haskey(get_settings(m), :pgap_value) #&& (haskey(m.settings, :gensys2) ? get_setting(m, :gensys2) : false)
        if get_setting(m, :pgap_type) ==:ngdp
            shocks, z0 = ngdp_forecast_init(m, shocks, z0, cond_type = cond_type)
        end
    end

    # Get variables necessary to enforce the zero lower bound in the forecast
    ind_r = m.observables[:obs_nominalrate]
    ind_r_sh = m.exogenous_shocks[:rm_sh]
    zlb_value = forecast_zlb_value(m)

    if isa(system, RegimeSwitchingSystem)
        forecast(m, system, z0, shocks; cond_type = cond_type, enforce_zlb = enforce_zlb,
                 ind_r = ind_r, ind_r_sh = ind_r_sh, zlb_value = zlb_value)
    else
        forecast(system, z0, shocks; enforce_zlb = enforce_zlb,
                 ind_r = ind_r, ind_r_sh = ind_r_sh, zlb_value = zlb_value)
    end
end

function forecast(system::System{S}, z0::Vector{S},
    shocks::Matrix{S}; enforce_zlb::Bool = false, ind_r::Int = -1,
    ind_r_sh::Int = -1, zlb_value::S = 0.13/4) where {S<:AbstractFloat}

    # Unpack system
    T, R, C = system[:TTT], system[:RRR], system[:CCC]
    Q, Z, D = system[:QQ], system[:ZZ], system[:DD]
    Z_pseudo, D_pseudo = system[:ZZ_pseudo], system[:DD_pseudo]

    # Setup
    nshocks = size(R, 2)
    nstates = size(T, 2)
    nobs    = size(Z, 1)
    npseudo = size(Z_pseudo, 1)
    horizon = size(shocks, 2)

    # Define our iteration function
    function iterate(z_t1, ϵ_t)
        z_t = C + T*z_t1 + R*ϵ_t

        # Change monetary policy shock to account for 0.13 interest rate bound
        if enforce_zlb
            interest_rate_forecast = getindex(D + Z*z_t, ind_r)
            if interest_rate_forecast < zlb_value
                # Solve for interest rate shock causing interest rate forecast to be exactly ZLB
                ϵ_t[ind_r_sh] = 0.
                z_t = C + T*z_t1 + R*ϵ_t
                ϵ_t[ind_r_sh] = getindex((zlb_value - D[ind_r] - Z[ind_r, :]'*z_t) / (Z[ind_r, :]' * R[:, ind_r_sh]), 1)

                # Forecast again with new shocks
                z_t = C + T*z_t1 + R*ϵ_t

                # Confirm procedure worked
                interest_rate_forecast = getindex(D + Z*z_t, ind_r)
                @assert interest_rate_forecast >= zlb_value - 0.01 "interest_rate_forecast = $interest_rate_forecast must be >= zlb_value - 0.01 = $(zlb_value - 0.01)"
            end
        end
        return z_t, ϵ_t
    end

    # Iterate state space forward
    states = zeros(S, nstates, horizon)
    states[:, 1], shocks[:, 1] = iterate(z0, shocks[:, 1])
    for t in 2:horizon
        states[:, t], shocks[:, t] = iterate(states[:, t-1], shocks[:, t])
    end

    # Apply measurement and pseudo-measurement equations
    obs    = D .+ Z*states
    pseudo = D_pseudo .+ Z_pseudo * states

    # Return forecasts
    return states, obs, pseudo, shocks
end

function forecast(m::AbstractDSGEModel, system::RegimeSwitchingSystem{S}, z0::Vector{S},
    shocks::Matrix{S}; cond_type::Symbol = :none, enforce_zlb::Bool = false, ind_r::Int = -1,
    ind_r_sh::Int = -1, zlb_value::S = 0.13/4) where {S<:AbstractFloat}

    # Determine how many regimes occur in the forecast, depending
    # on whether we need to subtract conditional forecast regimes or not
    n_fcast_reg = get_setting(m, :n_fcast_regimes)
    if cond_type != :none
        n_fcast_reg -= get_setting(m, :reg_post_conditional_end) - # Remove number of regimes switched since starting
            get_setting(m, :reg_forecast_start) - 1                # the conditional forecast
        if n_fcast_reg == 0 # Then no new regimes after conditional forecast ends
            n_fcast_reg = 1 # So we use the last conditional forecast regime
        end
    end

    Ts = Vector{Matrix{Float64}}(undef, n_fcast_reg)
    Rs = Vector{Matrix{Float64}}(undef, n_fcast_reg)
    Cs = Vector{Vector{Float64}}(undef, n_fcast_reg)
    Qs = Vector{Matrix{Float64}}(undef, n_fcast_reg)
    Zs = Vector{Matrix{Float64}}(undef, n_fcast_reg)
    Ds = Vector{Vector{Float64}}(undef, n_fcast_reg)
    Z_pseudos = Vector{Matrix{Float64}}(undef, n_fcast_reg)
    D_pseudos = Vector{Vector{Float64}}(undef, n_fcast_reg)

    # Determine in which regime the forecast starts, after accounting for
    # conditional forecast regimes, if applicable
    reg_fcast_cond_start = (cond_type == :none) ? get_setting(m, :reg_forecast_start) :
        max(get_setting(m, :reg_forecast_start), get_setting(m, :reg_post_conditional_end))

    # Unpack system
    for (ss_ind, sys_ind) in enumerate(reg_fcast_cond_start:get_setting(m, :n_regimes))
        Ts[ss_ind], Rs[ss_ind], Cs[ss_ind] = system[sys_ind, :TTT], system[sys_ind, :RRR], system[sys_ind, :CCC]
        Qs[ss_ind], Zs[ss_ind], Ds[ss_ind] = system[sys_ind, :QQ], system[sys_ind, :ZZ], system[sys_ind, :DD]
        Z_pseudos[ss_ind], D_pseudos[ss_ind] = system[sys_ind, :ZZ_pseudo], system[sys_ind, :DD_pseudo]
    end

    # Setup
    nshocks = size(Rs[1], 2)
    nstates = size(Ts[1], 2)
    nobs    = size(Zs[1], 1)
    npseudo = size(Z_pseudos[1], 1)
    horizon = size(shocks, 2)

    # Define our iteration function
    function iterate(z_t1, ϵ_t, T, R, C, Q, Z, D)
        z_t = C + T*z_t1 + R*ϵ_t

        # Change monetary policy shock to account for 0.13 interest rate bound
        if enforce_zlb
            interest_rate_forecast = getindex(D + Z*z_t, ind_r)
            if interest_rate_forecast < zlb_value
                # Solve for interest rate shock causing interest rate forecast to be exactly ZLB
                ϵ_t[ind_r_sh] = 0. # get forecast when MP shock
                z_t = C + T*z_t1 + R*ϵ_t # is zeroed out
                z_t_old = C + T*z_t1 + R*ϵ_t
                ϵ_t[ind_r_sh] = getindex((zlb_value - D[ind_r] - Z[ind_r, :]'*z_t) / (Z[ind_r, :]' * R[:, ind_r_sh]), 1)

                # Forecast again with new shocks
                z_t = C + T*z_t1 + R*ϵ_t

                # Confirm procedure worked
                interest_rate_forecast = getindex(D + Z*z_t, ind_r)

                # Subtract a small number to deal with numerical imprecision
                @assert interest_rate_forecast >= zlb_value - 0.01 "interest_rate_forecast = $interest_rate_forecast must be >= zlb_value - 0.01 = $(zlb_value - 0.01)."
            end
        end
        return z_t, ϵ_t
    end

    # These are the fcast regime indices
    regime_inds = get_fcast_regime_inds(m, horizon, cond_type)

    # Iterate state space forward
    states = zeros(nstates, horizon)
    obs = zeros(nobs, horizon)
    pseudo = zeros(npseudo, horizon)

    states[:, 1], shocks[:, 1] = iterate(z0, shocks[:, 1], Ts[1], Rs[1], Cs[1], Qs[1], Zs[1], Ds[1])
    obs[:, 1] = Ds[1] .+ Zs[1]*states[:, 1]
    pseudo[:, 1] = D_pseudos[1] .+ Z_pseudos[1] * states[:, 1]

    # If there's multiple regimes in forecast period, go through each set of indices. Otherwise, just take the first set
    if length(regime_inds) > 1
        for i in 1:length(regime_inds)
            ts = regime_inds[i]
            if ts[1] == 1 # We already did this period in the initialization
                ts = 2:ts[end]
            end
            for t in ts
                states[:, t], shocks[:, t] = iterate(states[:, t - 1], shocks[:, t], Ts[i], Rs[i], Cs[i], Qs[i], Zs[i], Ds[i])
                obs[:, t] = Ds[i] .+ Zs[i] * states[:, t]
                pseudo[:, t] = D_pseudos[i] .+ Z_pseudos[i] * states[:, t]
            end
        end
    else
        for t in 2:horizon
            states[:, t], shocks[:, t] = iterate(states[:, t - 1], shocks[:, t], Ts[1], Rs[1], Cs[1], Qs[1], Zs[1], Ds[1])
            obs[:, t] = Ds[1] .+ Zs[1] * states[:, t]
            pseudo[:, t] = D_pseudos[1] .+ Z_pseudos[1] * states[:, t]
        end
    end

    # Return forecasts
    return states, obs, pseudo, shocks
end

function get_fcast_regime_inds(m::AbstractDSGEModel, horizon::Int, cond_type::Symbol)

    fcast_start_date = (cond_type == :none) ? date_forecast_start(m) : # If conditional forecast, then we want to start forecasting
        max(date_forecast_start(m), iterate_quarters(date_conditional_end(m), 1)) # from the period after the conditional end period
    last_date = iterate_quarters(fcast_start_date, -1) # minus 1 to get one date before the forecast start
    last_ind = 1

    # Construct vector of time periods for each regime in the forecast periods after the conditional forecast
    regime_inds = Vector{UnitRange{Int}}(undef, 0)
    for i in 1:(get_setting(m, :n_regimes) - 1) # -1 because, we add the last regime:horizon to end
        if get_setting(m, :regime_dates)[i] < fcast_start_date # date_forecast_start(m)
            continue
        else
            qtr_diff = subtract_quarters(get_setting(m, :regime_dates)[i], last_date)
            regime_inds = push!(regime_inds, last_ind:(last_ind + qtr_diff - 1))
            last_ind = last_ind + qtr_diff
            last_date = get_setting(m, :regime_dates)[i]
        end
    end
    regime_inds = push!(regime_inds, last_ind:horizon) # This covers case where final history regime is also first (and only) forecast regime.
end
