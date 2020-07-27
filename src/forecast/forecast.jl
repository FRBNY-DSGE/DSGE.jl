"""
```
forecast(m, system, z0; enforce_zlb = false, shocks = Matrix{S}(undef, 0,0))

forecast(m, altpolicy, z0, shocks)

forecast(system, z0, shocks; enforce_zlb = false)
```

The first method produces a forecast, given a state space system,
initial state, and shocks, using information about the desired forecast
contained in `m`. It enforces the ZLB by using monetary policy shocks.

The second method is similar but differs in two ways. First,
it produces forecasts specifically when an alternative policy
is used. Second, it enforces the ZLB by treating it as a temporary
alternative policy.

The third method is an internal function used by the first two methods.

### Inputs

- `system::System{S}`: state-space system matrices
- `z0::Vector{S}`: state vector in the final historical period
- `shocks::Matrix{S}`: `nshocks` x `nperiods` matrix of shocks to use when
  forecasting. Note that in the first method, `nperiods` doesn't necessarily
  have to equal `forecast_horizons(m)`; it will be truncated or padded with
  zeros appropriately

**Method 1 and 2 only:**

- `m::AbstractDSGEModel`

**Method 2 only:**

- `altpolicy::Symbol`: Which alternative policy is being used
- `obs::Matrix{S <: Real}`: matrix of forecasted observables

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

**Method 2 only:**

- `temporary_altpolicy_max_iter::Int`: maximum number of iterations for which the
    function attempts to enforce the ZLB as a temporary alternative policy. Defaults to 10.

- `set_zlb_regime_vals::Function`: user-provided function that adds additional regimes to
    regime-switching parameters if not enough regimes exist to impose the ZLB
    as a temporary alternative policy. Defaults to `identity`, and nothing will happen
    if this is the case.

- `tol::{<: Real}`: Tolerance for the smallest permissible value for the nominal interest rate.
    Defaults to -1e-14.

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
            # If Regime Switching System, need to get proper QQ in each forecast period.
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

    # TODO: check that these handling of the altpolicy cases work properly
    if (alt_policy.solve != identity &&
        alt_policy.forecast_init != identity)
        shocks, z0 = alt_policy.forecast_init(m, shocks, z0, cond_type = cond_type)
    end

#=
    # This code block has been moved to forecast_one_draw in drivers.jl
    # b/c we want to actually alter the smoothed states so that pgap or ygap
    # starts from the user-specified gap value and that this start is actually recorded.
    # TODO: check that these handling of the altpolicy cases work properly
    # Added separately to cover case where you're not using altpolicy but are using gensys2
    if haskey(m.settings, :pgap_type) && haskey(get_settings(m), :pgap_value)
        if get_setting(m, :pgap_type) == :ngdp
            _, z0 = ngdp_forecast_init(m, shocks, z0, cond_type = cond_type)
        end
    end

    if haskey(m.settings, :ygap_type) && haskey(get_settings(m), :ygap_value)
        if get_setting(m, :ygap_type) == :smooth_ait_gdp
            _, z0 = smooth_ait_gdp_forecast_init(m, shocks, z0, cond_type = cond_type)
        end
    end
=#

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
            get_setting(m, :reg_forecast_start)                    # the conditional forecast
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

    # n_regimes is incorrectly being computed, as is n_fcast_reg
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
    obs[:, 1] = Ds[1] .+ Zs[1] * states[:, 1]
    pseudo[:, 1] = D_pseudos[1] .+ Z_pseudos[1] * states[:, 1]

    # If there's multiple regimes in forecast period, go through each set of indices. Otherwise, just take the first set
    if length(regime_inds) > 1
        for i in 1:length(regime_inds)
            ts = regime_inds[i]
            if length(ts) > 0 # Check in case the regime is empty
                if ts[1] == 1 # We already did this period in the initialization
                    ts = 2:ts[end]
                end
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

# Automatic enforcing of the ZLB as a temporary alternative policy
function forecast(m::AbstractDSGEModel, altpolicy::Symbol, z0::Vector{S}, states::AbstractMatrix{S},
                  obs::AbstractMatrix{S}, pseudo::AbstractMatrix{S}, shocks::AbstractMatrix{S};
                  cond_type::Symbol = :none, uncertain_altpolicy::Bool = false,
                  temporary_altpolicy_max_iter::Int = 10, set_zlb_regime_vals::Function = identity,
                  tol::S = -1e-14) where {S <: Real}

    # Grab "original" settings" so they can be restored later
    is_regime_switch = haskey(get_settings(m), :regime_switching) ? get_setting(m, :regime_switching) : false
    is_replace_eqcond = haskey(get_settings(m), :replace_eqcond) ? get_setting(m, :replace_eqcond) : false
    is_gensys2 = haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false

    # Grab some information about the forecast
    n_hist_regimes = haskey(DSGE.get_settings(m), :n_hist_regimes) ? get_setting(m, :n_hist_regimes) : 1
    if haskey(DSGE.get_settings(m), :regime_dates) # dates of first and last regimes we need to add for the alt policy
        if length(get_setting(m, :regime_dates)) > 1
            start_altpol_date = get_setting(m, :regime_dates)[2]
            end_altpol_date   = get_setting(m, :regime_dates)[get_setting(m, :reg_post_conditional_end)] +
                Dates.Month(3) * (forecast_horizons(m; cond_type = cond_type) + 1)
        else
            start_altpol_date = date_forecast_start(m)
            end_altpol_date   = iterate_quarters(start_altpol_date, forecast_horizons(m; cond_type = cond_type) + 1)
        end
    else
        start_altpol_date = date_forecast_start(m)
        end_altpol_date   = iterate_quarters(start_altpol_date, forecast_horizons(m; cond_type = cond_type) + 1)
    end

    # Figure out which periods need temporary ZLB regimes
#     which_zlb_regimes = findall(obs[get_observables(m)[:obs_nominalrate], :] .<
#                                 get_setting(m, :forecast_zlb_value)) .+ n_hist_regimes # add historical regimes to get the right regime number
#     if is_regime_switch && cond_type != :none
#         which_zlb_regimes .+= get_setting(m, :reg_post_conditional_end) - # additional conditional regimes should be added
#             get_setting(m, :reg_forecast_start)
#     end

    first_zlb_regime = findfirst(obs[get_observables(m)[:obs_nominalrate], :] .<
                                 get_setting(m, :forecast_zlb_value))
    if !isnothing(first_zlb_regime) # Then there are ZLB regimes to enforce
        first_zlb_regime += n_hist_regimes
        if is_regime_switch && cond_type != :none
            first_zlb_regime += get_setting(m, :reg_post_conditional_end) - get_setting(m, :reg_forecast_start)
        end

        # Setup for loop enforcing zero rate
        states = similar(states)
        pseudo = similar(pseudo)
        obs = copy(obs) # This way, we don't overwrite the underlying obs matrix

        # Now iteratively impose ZLB periods until there are no negative rates in the forecast horizon
        for iter in 0:(size(obs, 2) - 1)
            # Calculate the number of ZLB regimes. For now, we add in a separate regime for every
            # period b/n the first and last ZLB regime in the forecast horizon. It is typically the case
            # that this is necessary anyway but not always, especially depending on the drawn shocks
            n_total_regimes = first_zlb_regime + iter + 1 # plus 1 for lift off

            # Set up parameters if there are switching parameter values.
            #
            # User needs to provide a function which takes in the model object `m`
            # and the total number of regimes (after adding the required temporary regimes),
            # and sets up regime-switching parameters for these new additional regimes.
            if set_zlb_regime_vals != identity
                set_zlb_regime_vals(m, n_total_regimes)
            end

            # Set up regime dates
            altpol_regime_dates = Dict{Int, Date}(1 => date_presample_start(m))
            for (regind, date) in zip(2:n_total_regimes, # We want to have enough regimes to fill up altpol_regime_dates properly
                                      start_altpol_date:Dates.Month(3):end_altpol_date)
                altpol_regime_dates[regind] = date
            end
            m <= Setting(:regime_dates, altpol_regime_dates)
            m <= Setting(:regime_switching, true)
            setup_regime_switching_inds!(m)

            # Set up replace_eqcond entries
            m <= Setting(:replace_eqcond, true)
            m <= Setting(:gensys2, true)
            replace_eqcond = Dict{Int, Function}()     # Which rule to replace with in which periods
            for regind in first_zlb_regime:(n_total_regimes - 1)
                replace_eqcond[regind] = DSGE.zero_rate_replace_eq_entries # Temp ZLB rule in this regimes
            end

            replace_eqcond[n_total_regimes] = if altpolicy == :ait
                ait_replace_eq_entries # Add permanent AIT regime
            elseif altpolicy == :ngdp
                ngdp_replace_eq_entries # Add permanent NGDP regime
            elseif altpolicy == :smooth_ait_gdp
                smooth_ait_gdp_replace_eq_entries # Add permanent smooth AIT-GDP regime
            elseif altpolicy == :smooth_ait_gdp_alt
                smooth_ait_gdp_alt_replace_eq_entries # Add permanent smooth AIT-GDP (alternative specification) regime
            else # Default to historical policy
                eqcond
            end
            m <= Setting(:replace_eqcond_func_dict, replace_eqcond)

            # Recompute to account for new regimes
            system = compute_system(m; apply_altpolicy = true)

            # Forecast!
            states[:, :], obs[:, :], pseudo[:, :] = forecast(m, system, z0; cond_type = cond_type, shocks = shocks)

            if all(obs[get_observables(m)[:obs_nominalrate], :] .> tol)
                break
            end
        end
    end

    if !all(obs[get_observables(m)[:obs_nominalrate], :] .> tol)
        @warn "Unable to enforce the ZLB."
    end

    # Restore original settings
    m <= Setting(:regime_switching, is_regime_switch)
    m <= Setting(:replace_eqcond,   is_replace_eqcond)
    m <= Setting(:gensys2,          is_gensys2)

    return states, obs, pseudo
end

function get_fcast_regime_inds(m::AbstractDSGEModel, horizon::Int, cond_type::Symbol)

    fcast_post_cond_date = (cond_type == :none) ? date_forecast_start(m) : # If cond forecast, then we want to start forecasting
        max(date_forecast_start(m), iterate_quarters(date_conditional_end(m), 1)) # from the period after the conditional end period
    # last_date = iterate_quarters(fcast_post_cond_date, -1) # minus 1 to get one date before the forecast start
    last_date = fcast_post_cond_date
    last_ind = 1

    # Construct vector of time periods for each regime in the forecast periods after the conditional forecast
    regime_inds = Vector{UnitRange{Int}}(undef, 0)
    start_reg = (cond_type == :none) ? get_setting(m, :reg_forecast_start) :
        max(get_setting(m, :reg_forecast_start), get_setting(m, :reg_post_conditional_end))
    for i in start_reg:(get_setting(m, :n_regimes) - 1) # Last regime handled separately
        qtr_diff = subtract_quarters(get_setting(m, :regime_dates)[i + 1], last_date)
        push!(regime_inds, last_ind:(last_ind + qtr_diff - 1))
        last_ind += qtr_diff
        last_date = get_setting(m, :regime_dates)[i + 1]
    end
    regime_inds = push!(regime_inds, last_ind:horizon) # Covers case where final hist regime is also first (and only) forecast regime.

    return regime_inds
end
