"""
```
forecast(m, system, z0; cond_type = :none,
    enforce_zlb = false, shocks = Matrix{S}(undef, 0,0))

forecast(m, altpolicy, z0, states, obs, pseudo, shocks)

forecast(system, z0, shocks; enforce_zlb = false)

forecast(m, system, z0, shocks; cond_type = :none, enforce_zlb = false)
```

The first method produces a forecast, given a state space system,
initial state, and shocks, using information about the desired forecast
contained in `m`. It enforces the ZLB by using monetary policy shocks.

The second method is similar but differs in two ways. First,
it produces forecasts specifically when an alternative policy
is used. Second, it enforces the ZLB by treating it as a temporary
alternative policy.

The third and fourth methods are internal functions used by the first two methods.

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
                ind_ant1 = m.exogenous_shocks[Symbol(get_setting(m, :monetary_policy_shock), :l1)]
                ind_antn = m.exogenous_shocks[Symbol(get_setting(m, :monetary_policy_shock), "l$(n_mon_anticipated_shocks(m))")]
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
#=    alt_policy = alternative_policy(m)

    # TODO: check that these handling of the altpolicy cases work properly
    if (alt_policy.solve != identity &&
        alt_policy.forecast_init != identity && (!haskey(m.settings, :initialize_pgap_ygap) || get_setting(m, :initialize_pgap_ygap)))
        shocks, z0 = alt_policy.forecast_init(m, shocks, z0, cond_type = cond_type)
    end=#

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
    ind_r = m.observables[get_setting(m, :nominal_rate_observable)]
    ind_r_sh = m.exogenous_shocks[get_setting(m, :monetary_policy_shock)]
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
    ind_r_sh::Int = -1, zlb_value::S = 0.1/4) where {S<:AbstractFloat}

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
    ind_r_sh::Int = -1, zlb_value::S = 0.1/4) where {S<:AbstractFloat}

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
    check_zero_rate = (haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false) &&
        (:zero_rate in map(x -> x.alternative_policy.key, collect(values(get_setting(m, :regime_eqcond_info)))))
    function iterate(z_t1, ϵ_t, T, R, C, Q, Z, D)
        z_t = C + T*z_t1 + R*ϵ_t

        # Change monetary policy shock to account for 0.13 interest rate bound
        if enforce_zlb
            interest_rate_forecast = getindex(D + Z*z_t, ind_r)
            if interest_rate_forecast < zlb_value
                continue_enforce = check_zero_rate ? abs.(Z[ind_r, :]' * R[:, ind_r_sh]) > 1e-4 : true # make sure we don't actually divide by zero (which may occur if R[:, ind_r_sh] .== 0.
                if continue_enforce
                    # Solve for interest rate shock causing interest rate forecast to be exactly ZLB
                    ϵ_t[ind_r_sh] = 0. # get forecast when MP shock
                    z_t = C + T*z_t1 + R*ϵ_t # is zeroed out
                    z_t_old = C + T*z_t1 + R*ϵ_t
                    ϵ_t[ind_r_sh] = getindex((zlb_value - D[ind_r] - Z[ind_r, :]'*z_t) / (Z[ind_r, :]' * R[:, ind_r_sh]), 1)

                    # Forecast again with new shocks
                    z_t = C + T*z_t1 + R*ϵ_t

                    # Confirm procedure worked
                    interest_rate_forecast = getindex(D + Z*z_t, ind_r)
                    if isnan(interest_rate_forecast)
                        ϵ_t[ind_r_sh] = 0. # get forecast when MP shock
                        z_t = C + T*z_t1 + R*ϵ_t # is zeroed out
                        z_t_old = C + T*z_t1 + R*ϵ_t
                        ϵ_t[ind_r_sh] = getindex((zlb_value - D[ind_r] - Z[ind_r, :]'*z_t) / (Z[ind_r, :]' * R[:, ind_r_sh]), 1)
                    end

                    # Subtract a small number to deal with numerical imprecision
                    @assert interest_rate_forecast >= zlb_value - 0.01 "interest_rate_forecast = $interest_rate_forecast must be >= zlb_value - 0.01 = $(zlb_value - 0.01)."
                end
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
function forecast(m::AbstractDSGEModel, z0::Vector{S}, states::AbstractMatrix{S},
                  obs::AbstractMatrix{S}, pseudo::AbstractMatrix{S}, shocks::AbstractMatrix{S};
                  cond_type::Symbol = :none, set_zlb_regime_vals::Function = identity,
                  set_info_sets_altpolicy::Function = auto_temp_altpolicy_info_set,
                  update_regime_eqcond_info!::Function = default_update_regime_eqcond_info!,
                  tol::S = -1e-5, rerun_smoother::Bool = false,
                  df = nothing, draw_states::Bool = false,
                  nan_failures::Bool = false,
                  histstates::AbstractMatrix{S} = Matrix{S}(undef, 0, 0),
                  histshocks::AbstractMatrix{S} = Matrix{S}(undef, 0, 0),
                  histpseudo::AbstractMatrix{S} = Matrix{S}(undef, 0, 0),
                  initial_states::AbstractVector{S} = Vector{S}(undef, 0)) where {S <: Real}

    # Grab "original" settings" so they can be restored later
    altpol = alternative_policy(m)
    is_regime_switch = haskey(get_settings(m), :regime_switching) ? get_setting(m, :regime_switching) : false
    is_replace_eqcond = haskey(get_settings(m), :replace_eqcond) ? get_setting(m, :replace_eqcond) : false
    is_gensys2 = haskey(get_settings(m), :gensys2) ? get_setting(m, :gensys2) : false
    original_info_set = haskey(get_settings(m), :tvis_information_set) ? get_setting(m, :tvis_information_set) : UnitRange{Int64}[]
    original_eqcond_dict = haskey(get_settings(m), :regime_eqcond_info) ? get_setting(m, :regime_eqcond_info) :
        Dict{Int, EqcondEntry}()

    # Information set - start of awareness of ZLB
    first_aware = if haskey(get_settings(m), :tvis_information_set)
        findfirst([maximum(original_info_set[i]) == get_setting(m, :n_regimes) for i in keys(original_info_set)])
    elseif haskey(get_settings(m), :regime_eqcond_info)
        findfirst([original_info_set[i].key == :zero_rate for i in keys(original_info_set)])
    else
        nothing
    end

    # Grab some information about the forecast
    n_hist_regimes = haskey(get_settings(m), :n_hist_regimes) ? get_setting(m, :n_hist_regimes) : 1
    has_reg_dates = haskey(get_settings(m), :regime_dates) # calculate dates of first & last regimes we need to add for the alt policy
    if has_reg_dates && length(get_setting(m, :regime_dates)) > 1
        start_altpol_reg  = get_setting(m, :reg_forecast_start)
        start_altpol_date = get_setting(m, :regime_dates)[start_altpol_reg]
        end_altpol_date   = date_forecast_end(m)
    else
        start_altpol_reg  = n_hist_regimes + 1
        start_altpol_date = date_forecast_start(m)
        end_altpol_date   = iterate_quarters(start_altpol_date, forecast_horizons(m; cond_type = cond_type) - 1)
    end
    altpol_reg_qtrrange   = quarter_range(start_altpol_date, end_altpol_date)
    n_altpol_reg_qtrrange = length(altpol_reg_qtrrange)

    # Start imposing ZLB at the quarter before liftoff quarter
    has_neg_rates    = view(obs, get_observables(m)[:obs_nominalrate], :) .<= forecast_zlb_value(m) / 4.
    last_zlb_regime  = findlast(has_neg_rates)

    if !isnothing(last_zlb_regime) # Then there are ZLB regimes to enforce

        endo_success = false

        # Check original_eqcond_dict if any regimes use the zero rate rule.
        for (reg, v) in original_eqcond_dict
            if v.alternative_policy.key == :zero_rate &&
                (reg >= (cond_type == :none ? get_setting(m, :reg_forecast_start) : get_setting(m, :reg_post_conditional_end)))
                @warn "Regime $reg of regime_eqcond_info used zero_rate--to avoid gensys errors in computing the endogenous zlb, this regime is now being set to use $(altpol.key)."
                v.alternative_policy = altpol
            end
        end

        # Update first_zlb_regime to actually be a regime number, rather than index of forecasted obs
        pre_fcast_regimes = n_hist_regimes
        if cond_type != :none
            if is_regime_switch
                pre_fcast_regimes += get_setting(m, :reg_post_conditional_end) - get_setting(m, :reg_forecast_start)
            else
                pre_fcast_regimes += subtract_quarters(get_setting(m, :date_conditional_end),
                                                      get_setting(m, :date_forecast_start)) + 1
            end
        end

        # TODO: check if this behavior is really what we want, rather than imposing the ZLB regime
        # in the first period for which negative rates appear, i.e.
        # first_zlb_regime = findfirst(has_neg_rates)
        first_zlb_regime = pre_fcast_regimes + 1 # we start the ZLB after the history/conditional horizon

        # Update first awareness period if nothing before
        if isnothing(first_aware)
            first_aware = first_zlb_regime
        end

        # The temporary_altpol_length is the no. of regimes
        # in the history/conditional horizon that are ZLB.
        # m <= Setting(:temporary_altpol_length, last_zlb_regime - first_zlb_regime + 1)

        # Setup for loop enforcing zero rate
        orig_regimes      = haskey(get_settings(m), :n_regimes) ? get_setting(m, :n_regimes) : 1
        orig_regime_dates = haskey(get_settings(m), :regime_dates) ? get_setting(m, :regime_dates) : Dict{Int, Date}()
        orig_temp_zlb     = haskey(get_settings(m), :temporary_altpol_length) ? get_setting(m, :temporary_altpol_length) : nothing

        # Now iteratively impose ZLB periods until there are no negative rates in the forecast horizon

        # Set up for binary search
        max_zlb_regimes = haskey(get_settings(m), :max_temporary_altpol_length) ?
            get_setting(m, :max_temporary_altpol_length) - 1 : size(obs, 2) - 3 # subtract 1 b/c will add 1 later (see line 459)
        max_zlb_regimes += pre_fcast_regimes # so max_zlb_regimes is the regime number
        iter = min(last_zlb_regime, max_zlb_regimes)
        to_return = false
        high = max_zlb_regimes
        low = pre_fcast_regimes

        ## Binary search for correct number of zlb
        while true
            # Calculate the number of ZLB regimes. For now, we add in a separate regime for every
            # period b/n the first and last ZLB regime in the forecast horizon. It is typically the case
            # that this is necessary anyway but not always, especially depending on the drawn shocks
            @show iter
            n_total_regimes = iter + 1 # plus 1 for lift off
            m <= Setting(:n_regimes, max(orig_regimes, n_total_regimes))

            ## TODO: Handle edge case where no. of ZLB needed is 0
	        m <= Setting(:temporary_altpol_length, iter - first_zlb_regime + 1) # 1 for first_zlb_regime, iter for each additional regime
            # Test if more/less ZLB Needed
            # Set up regime dates
            altpol_regime_dates = Dict{Int, Date}(1 => date_presample_start(m))
            if is_regime_switch # Add historical regimes
                for regind in 2:n_hist_regimes
                    altpol_regime_dates[regind] = orig_regime_dates[regind]
                end
            end

            # Ensure we don't accidentally tell the code to forecast beyond the forecast horizon
            altpol_reg_range = start_altpol_reg:n_total_regimes
            if length(altpol_reg_range) > n_altpol_reg_qtrrange
                if all(obs[get_observables(m)[:obs_nominalrate], :] .>
                       forecast_zlb_value(m) / 4.0 + tol) # Ensure we do not accidentally
                    obs[get_observables(m)[:obs_nominalrate], 1] = tol - 2. # assume success at enforcing ZLB (should be negative)
                end
                break
            end
            for (regind, date) in zip(altpol_reg_range, altpol_reg_qtrrange)
                altpol_regime_dates[regind] = date
            end

            if haskey(get_settings(m), :cred_vary_until) && isempty(altpol_reg_range) ? true : get_setting(m, :cred_vary_until) >= maximum(altpol_reg_range)
                m <= Setting(:regime_switching, true)
            else
                m <= Setting(:regime_dates, altpol_regime_dates)
                m <= Setting(:regime_switching, true)
                setup_regime_switching_inds!(m; cond_type = cond_type)
            end

            # Set up replace_eqcond entries
            m <= Setting(:replace_eqcond, true)
            m <= Setting(:gensys2, true)

            # In general, the update_regime_eqcond_info! function should not change the
            # EqcondEntry during historical/conditional horizon regimes, but any other
            # regimes in the forecast horizon should be/can be set.
            #
            # From above, we also assume that no ZLBs occur in the eqcond_dict that
            # is passed to update_regime_eqcond_info! and assume that
            # if there is a particular sequence of alternative policies the user
            # wants to occur after the ZLB, then the user will have incorporated
            # code that implements that sequence in their own update_regime_eqcond_info!
            # function. In the default DSGE policy, the regimes after the ZLB ends
            # are updated only if there is time-varying credibility
            # (specified by the Setting :cred_vary_until).
            update_regime_eqcond_info!(m, deepcopy(original_eqcond_dict), first_zlb_regime, n_total_regimes, altpol)

            # Set up parameters if there are switching parameter values.
            #
            # User needs to provide a function which takes in the model object `m`
            # and the total number of regimes (after adding the required temporary regimes),
            # and sets up regime-switching parameters for these new additional regimes.
#=            if set_zlb_regime_vals != identity &&
                (!haskey(get_settings(m), :cred_vary_until) || n_total_regimes > get_setting(m, :cred_vary_until))
                set_zlb_regime_vals(m, n_total_regimes)
            end=#
            if set_zlb_regime_vals != identity
                set_zlb_regime_vals(m, n_total_regimes)
            end

            # set up the information sets TODO: add checkfor whether or not we even need to update the tvis_info_set
            set_info_sets_altpolicy(m, get_setting(m, :n_regimes), first_aware)
            #=if haskey(get_settings(m), :cred_vary_until) && get_setting(m, :cred_vary_until) >= n_total_regimes
                set_info_sets_altpolicy(m, get_setting(m, :cred_vary_until) + 1, zlb_at_first)
            else
                set_info_sets_altpolicy(m, get_setting(m, :n_regimes), first_zlb_regime)
            end=#

            # TODO: maybe delete tvis, since tvis should always be on given the fact we always run set_info_sets
            # tvis = haskey(get_settings(m), :tvis_information_set) && !isempty(get_setting(m, :tvis_information_set))

            # Recompute to account for new regimes
            system = compute_system(m; tvis = true)

            if rerun_smoother # if state space system changes, then the smoothed states will also change generally
                histstates, histshocks, histpseudo, initial_states =
                    smooth(m, df, system; cond_type = cond_type, draw_states = draw_states)
                z0 = histstates[:, end]
            end

            # Forecast! Note that since forecast(m, system, z0) allocates new matrices for states, obs, and pseudo,
            # this step merely changes to which matrices the local variables states, obs, and pseudo point
            # w/in this function's closure. Thus, the matrices states, obs, and pseudo which are first
            # passed into this function will not be over-written.
            states, obs, pseudo = forecast(m, system, z0; cond_type = cond_type, shocks = shocks)

            # Successful endogenous bounding?
            endo_success = all(obs[get_observables(m)[:obs_nominalrate], :] .> get_setting(m, :zero_rate_zlb_value) / 4. + tol)

            if !endo_success && iter == max_zlb_regimes
                states, obs, pseudo = forecast(m, system, z0; cond_type = cond_type, shocks = shocks,
                                               enforce_zlb = true)
                endo_success = true
                break
            end

            # Delete extra regimes added to implement the temporary alternative policy, or else updating the parameters
            # in forecast_one will not work.
            if set_zlb_regime_vals != identity
                for p in m.parameters
                    if haskey(p.regimes, :value)
                        if length(p.regimes[:value]) > orig_regimes
                            for i in (orig_regimes + 1):length(p.regimes[:value])
                                delete!(p.regimes[:value], i)
                            end
                        end
                    end
                end
            end

            if to_return || (endo_success &&
                             (last_zlb_regime - 1 <= iter <= last_zlb_regime + 2)) ## We ran this iteration to return the answer.
                @assert endo_success "Binary search for endo ZLB is breaking even though endo_success is false"

                break
            end

            ## Note in below 2 is really :endogenous_zlb_lookback and 3 is :endogenous_zlb_lookahead
            ## Reset temporary_altpol_length in each case
            ## Return here refers to "break"
            ## TODO: Conditional period included here b/c we NaN out conditional nominal rates
            ### in setup_flexait_tempzlb! This contradicts the earlier setup where we start
            ### at reg_post_conditional_end instead of reg_forecast_start.

            # If more ZLB necessary
            @show iter, low, high, first_zlb_regime, last_zlb_regime, pre_fcast_regimes
            if !endo_success
#=                low = iter + 1
                iter_old = iter
                if iter == first_zlb_regime
                    iter = min(iter+3, max_zlb_regimes)=#
                low = iter + 1 # Set low to iter + 1 if endo_success fails
                if iter == last_zlb_regime
                    # Our initial guess of last_zlb_regime failed.
                    # We will try guessing 3 periods ahead (assuming we don't exceed max_zlb_regimes)
                    iter = min(iter + 3, max_zlb_regimes)
                elseif iter == max_zlb_regimes || low == high ## Second or is redundant
                    # We have hit the max number of allowed regimes for ZLB, so
                    # run Fixed ZLB with unanticipated shocks to enforce ZLB
                    break
                elseif iter == last_zlb_regime + 3
                    # Guessing 3 periods after last_zlb_regime failed,
                    # so let's try max_zlb_regimes
                    iter = max_zlb_regimes
                elseif iter == last_zlb_regime + 1
                    # Increment by one period b/c while
                    # last_zlb_regime doesn't enforce, max_zlb_regime or
                    # last_zlb_regime + 3 definitely do
                    iter = last_zlb_regime + 2
                elseif iter == last_zlb_regime + 2
                    # Increment by one period b/c while
                    # last_zlb_regime doesn't enforce,
                    # max_zlb_regime definitely does
                    iter = last_zlb_regime + 3
                    to_return = true # Return b/c last_zlb_regime + 3 will work
                elseif iter == last_zlb_regime - 2
                    # Increment by one period b/c while
                    # low = pre_fcast_regimes doesn't work
                    # but last_zlb_regime does
                    iter = last_zlb_regime - 1
                elseif iter == last_zlb_regime - 1
                    iter = last_zlb_regime
                    to_return = true # Return b/c last_zlb_regime will work
                else
                    # Continue Binary Search
                    new_iter = floor(Int64, (low + high) / 2)
                    if iter == new_iter
                        # Set to high b/c high must be working
                        # but low isn't (e.g high = low + 1)
                        iter = high
                        to_return = true
                        break
                    else
                        iter = new_iter
                    end
                end
            else
                # If ZLB works (no negative rates w/out using unanticipated mp shocks)
                high = iter # Set high to iter so we search for ZLB regimes less than or equal to iter
                if iter == last_zlb_regime
                    # Let's try a final ZLB regime 2 before last_zlb_regime, but above first_zlb_regime as bound
                    iter = max(last_zlb_regime - 2, first_zlb_regime)
                elseif iter == last_zlb_regime - 2
                    # If last_zlb_regime - 2 works, let's try the first possible ZLB regime value
                    iter = first_zlb_regime
                elseif iter == first_zlb_regime
                    # Minimal number of ZLB works, so we break here
                    # and use iter = first_zlb_regime
                    break
                elseif iter == last_zlb_regime - 1
                    # If this is true, we must have tried last_zlb_regime successfully,
                    # but values below or equal to last_zlb_regime - 2 failed.
                    # We will use last_zlb_regime - 1
                    break
                elseif iter == last_zlb_regime + 1
                    # last_zlb_regime had failed, but last_zlb_regime + 3 worked,
                    # so we use last_zlb_regime + 1
                    break
                elseif iter == last_zlb_regime + 2
                    # last_zlb_regime had failed, but last_zlb_regime + 3 worked,
                    # so we use last_zlb_regime + 2
                    break
                elseif iter == last_zlb_regime + 3
                    # last_zlb_regime had failed, but max_zlb_regimes + pre_fcast_regimes and this worked
                    # so we will guess last_zlb_regime + 1 to check if it works
                    iter = last_zlb_regime + 1
                elseif high == low || iter == low
                    # Return this
                    break
                elseif iter == max_zlb_regimes
                    # Continue Binary Search on [last_zlb_regime + 4, max_zlb_regimes]
                    iter = floor(Int64, (max_zlb_regimes + last_zlb_regime + 4) / 2)
                    high = max_zlb_regimes
                    low = last_zlb_regime + 4
                else
                    # Continue Binary Search
                    new_iter = floor(Int64, (low + high) / 2)
                    if iter == new_iter
                        # Set to high b/c high must be working
                        # but low isn't (e.g high = low + 1)
                        iter = high
                        to_return = true
                        break
                    else
                        iter = new_iter
                    end
                end
            end
        end
    else
        endo_success = true

        if rerun_smoother
            return states, obs, pseudo, histstates, histshocks, histpseudo, initial_states
        else
            return states, obs, pseudo
        end
    end

    if nan_failures && !endo_success
        @warn "Unable to enforce the ZLB. Throwing NaN for forecasts"
        states .= NaN
        obs    .= NaN
        pseudo .= NaN
    end

    # Restore original settings
    m <= Setting(:regime_switching, is_regime_switch)
    m <= Setting(:regime_dates, orig_regime_dates)
    if !isempty(orig_regime_dates)
        setup_regime_switching_inds!(m; cond_type = cond_type)
    end
    m <= Setting(:replace_eqcond,   is_replace_eqcond)
    m <= Setting(:gensys2,          is_gensys2)
    if !is_replace_eqcond
        delete!(get_settings(m), :regime_eqcond_info)
    end
    if isempty(original_eqcond_dict)
        delete!(get_settings(m), :regime_eqcond_info)
    else
        m <= Setting(:regime_eqcond_info, original_eqcond_dict)
    end
    if haskey(get_settings(m), :tvis_information_set)
        m <= Setting(:tvis_information_set, original_info_set)
    end
    if isnothing(orig_temp_zlb)
        delete!(get_settings(m), :temporary_altpol_length)
    else
        m <= Setting(:temporary_altpol_length, orig_temp_zlb)
    end

    if rerun_smoother
        # Need to return histstates, etc., to update the histstates, etc.,
        # which are returned in the forecast_output dictionary
        return states, obs, pseudo, histstates, histshocks, histpseudo, initial_states
    else
        return states, obs, pseudo
    end
end

function get_fcast_regime_inds(m::AbstractDSGEModel, horizon::Int, cond_type::Symbol;
                               start_index::Int = 0)

    fcast_post_cond_date = (cond_type == :none) ? date_forecast_start(m) : # If cond forecast, then we want to start forecasting
        max(date_forecast_start(m), iterate_quarters(date_conditional_end(m), 1)) # from the period after the conditional end period
    # last_date = iterate_quarters(fcast_post_cond_date, -1) # minus 1 to get one date before the forecast start
    last_date = fcast_post_cond_date
    last_ind = 1

    # Construct vector of time periods for each regime in the forecast periods after the conditional forecast
    start_reg   = (cond_type == :none) ? get_setting(m, :reg_forecast_start) :
        max(get_setting(m, :reg_forecast_start), get_setting(m, :reg_post_conditional_end))
    regime_inds = Vector{UnitRange{Int}}(undef, get_setting(m, :n_regimes) - start_reg + 1)
    for (rgi, i) in enumerate(start_reg:(get_setting(m, :n_regimes) - 1)) # Last regime handled separately
        qtr_diff = subtract_quarters(get_setting(m, :regime_dates)[i + 1], last_date)
        regime_inds[rgi] = (start_index + last_ind):(start_index + last_ind + qtr_diff - 1)
        last_ind += qtr_diff
        last_date = get_setting(m, :regime_dates)[i + 1]
    end
    regime_inds[end] = (start_index + last_ind):(start_index + horizon) # Covers case where final hist regime is also first (and only) forecast regime.

    return regime_inds
end

function auto_temp_altpolicy_info_set(m::AbstractDSGEModel, n_regimes::Int, first_zlb_regime::Int)
    tvis_infoset = Vector{UnitRange{Int}}(undef, n_regimes)
    for r in 1:n_regimes
        if r < first_zlb_regime
            tvis_infoset[r] = r:r
        else
            tvis_infoset[r] = r:n_regimes
        end
    end
    m <= Setting(:tvis_information_set, tvis_infoset)
end

# eqcond_dict should be a dictionary whose historical/conditional horizon regimes
# should already be set (and won't be affected by this function.
# The regimes which will be altered are those for which the temporary ZLB will apply
# (according to `zlb_start_regime` and `liftoff_regime`), and if there is
# time-varying credibility, forecast regimes past the ZLB.
#
# In general, the update_regime_eqcond_info! function should not change the
# EqcondEntry during historical/conditional horizon regimes, but any other
# regimes in the forecast horizon should be/can be set.
function default_update_regime_eqcond_info!(m::AbstractDSGEModel, eqcond_dict::AbstractDict{Int64, EqcondEntry},
                                            zlb_start_regime::Int64, liftoff_regime::Int64, liftoff_policy::AltPolicy)

    # TODO: add a check for making sure weights are properly set (e.g. enough regimes with weights to avoid error)
    # Populate eqcond_dict for regimes in which the ZLB should apply and update the weights as appropriate
    for regind in zlb_start_regime:(liftoff_regime - 1)
        if haskey(eqcond_dict, regind)
            # If regime index exists already, just directly change the
            # alternative_policy (which leaves the credibility weights unchanged)
            eqcond_dict[regind].alternative_policy = zero_rate() # Temp ZLB rule in this regimes
        else
            # If the regime index doesn't already exist, then add an EqcondEntry
            # for a ZLB that is perfectly credible
            weights = zeros(haskey(get_settings(m), :alternative_policies) ? length(get_setting(m, :alternative_policies)) + 1 : 1)
            weights[1] = 1.
            eqcond_dict[regind] = EqcondEntry(zero_rate(), weights)
        end
    end

    # Set the liftoff regime's EqcondEntry, using same logic as above
    if haskey(eqcond_dict, liftoff_regime)
        eqcond_dict[liftoff_regime].alternative_policy = liftoff_policy
    else
        weights = zeros(haskey(m.settings, :alternative_policies) ? length(get_setting(m, :alternative_policies)) + 1 : 1)
        weights[1] = 1.
        eqcond_dict[liftoff_regime] = EqcondEntry(liftoff_policy, weights)
    end

    # NOTE: the following block assumes there is only one temporary altpol
    # Update eqcond_dict to use the lift-off policy for any regimes between (inclusive)
    # the liftoff regime and the first regime for which credibility stops varying.
    if haskey(get_settings(m), :cred_vary_until) && get_setting(m, :cred_vary_until) >= liftoff_regime
        for z in liftoff_regime:(get_setting(m, :cred_vary_until) + 1)
            eqcond_dict[z].alternative_policy = liftoff_policy
        end
    end

    m <= Setting(:regime_eqcond_info, eqcond_dict)

    return m
end
