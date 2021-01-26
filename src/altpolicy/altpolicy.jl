"""
```
mutable struct AltPolicy
```

Type defining an alternative policy rule.

### Fields

- `key::Symbol`: alternative policy identifier

- `eqcond::Function`: a version of `DSGE.eqcond` which computes the equilibrium
  condition matrices under the alternative policy. Like `DSGE.eqcond`, it should
  take in one argument of mutable struct `AbstractDSGEModel` and return the `Γ0`,
  `Γ1`, `C`, `Ψ`, and `Π` matrices.

- `solve::Function`: a version of `DSGE.solve` which solves the model under the
  alternative policy. Like `DSGE.solve`, it should take in one argument of mutable
  struct `AbstractDSGEModel` and return the `TTT`, `RRR`, and `CCC` matrices.

- `forecast_init::Function`: a function that initializes forecasts under the
  alternative policy rule. Specifically, it accepts a model, an `nshocks` x
  `n_forecast_periods` matrix of shocks to be applied in the forecast, and a
  vector of initial states for the forecast. It must return a new matrix of
  shocks and a new initial state vector. If no adjustments to shocks or initial
  state vectors are necessary under the policy rule, this field may be omitted.

- `color::Colorant`: color to plot this alternative policy in. Defaults to blue.

- `linestyle::Symbol`: line style for forecast plots under this alternative
  policy. See options from `Plots.jl`. Defaults to `:solid`.

"""
mutable struct AltPolicy
    key::Symbol
    eqcond::Function
    solve::Function
    forecast_init::Function
    color::Colorant
    linestyle::Symbol
end

function AltPolicy(key::Symbol, eqcond_fcn::Function, solve_fcn::Function;
                   forecast_init::Function = identity,
                   color::Colorant = RGB(0., 0., 1.),
                   linestyle::Symbol = :solid)

    AltPolicy(key, eqcond_fcn, solve_fcn, forecast_init, color, linestyle)
end

Base.string(a::AltPolicy) = string(a.key)

function altpolicy_replace_eq_entries(key::Symbol)
    if key == :flexible_ait
        return flexible_ait_replace_eq_entries
    elseif key == :zero_rate
        return zero_rate_replace_eq_entries
    elseif key == :smooth_ait_gdp_alt
        smooth_ait_gdp_alt_replace_eq_entries
    elseif key == :smooth_ait_gdp
        smooth_ait_gdp_replace_eq_entries
    elseif key == :ait
        ait_replace_eq_entries
    elseif key == :ngdp
        ngdp_replace_eq_entries
    elseif key == :rw_zero_rate
        rw_zero_rate_replace_eq_entries
    elseif key == :rw
        rw_replace_eq_entries
    else
        error("AltPolicy $(key) has not been added to `altpolicy_replace_eq_entries` as having a method for replacing equilibrium conditions")
    end
end

"""
```
mutable struct EqcondEntry
```

Type to hold the entries in the regime_eqcond_info dictionary for
alternative policies, regime switching, and imperfect awareness.
"""
mutable struct EqcondEntry
    alternative_policy::Union{AltPolicy, Missing}
    weights::Union{Array{Float64, 1}, Missing}
end

function EqcondEntry()
    return EqcondEntry(missing, missing)
end

function EqcondEntry(altpol::AltPolicy)
    return EqcondEntry(altpol, missing)
end

"""
```
setup_permanent_altpol!(m::AbstractDSGEModel, altpol::AltPolicy;
    cond_type::Symbol = :none, remove_extra_eqcond_entry::Bool = true)
```

sets up `m` to implement a permanent alternative policy via exogenous regime-switching.
This function assumes that the settings `:date_forecast_start` and/or
`:date_conditional_end` are properly set. The latter setting only matters
if the keyword `cond_type` is set to `:semi` or `:full`.

The keyword `remove_extra_eqcond_entry` removes any `EqcondEntry`
entries in the setting `:regime_eqcond_info` after the
regime implementing the alternative policy, i.e. regimes in the forecast horizon.
If `false`, we will perform this step.
"""
function setup_permanent_altpol!(m::AbstractDSGEModel, altpol::AltPolicy;
                                 cond_type::Symbol = :none, remove_extra_eqcond_entry::Bool = true)

    # Figure out the date of the regime switch to the alt policy
    altpol_date = cond_type == :none ? date_forecast_start(m) : iterate_quarters(date_conditional_end(m), 1)

    # Add regime-switching setting
    if !haskey(get_settings(m), :regime_switching) ||
        (haskey(get_settings(m), :regime_switching) && !get_setting(m, :regime_switching))
        @warn "The setting :regime_switching is either missing or false. Changing it to true"
        m <= Setting(:regime_switching, true)
    end

    # Update regime dates
    if haskey(get_settings(m), :regime_dates)
        n_regs      = length(get_setting(m, :regime_dates)) # calculate length here rather than n_regimes in case
                                                            # setup_regime_switching_inds! has not been called already
        counter     = 0
        reg_missing = true                                  # if false, then altpol date is already in regime_dates
        reg_dates   = get_setting(m, :regime_dates)
        for reg in 1:n_regs
            if reg_dates[reg] > altpol_date
                break
            elseif reg_dates[reg] == altpol_date
                reg_missing = false
                break
            end
            counter += 1
        end
        if counter == n_regs # Case where all regimes occur in history
            altpol_reg = n_regs + 1
            reg_dates[altpol_reg] = altpol_date
        elseif reg_missing
            altpol_reg = counter + 1 # add 1 to get the regime for the alternative policy
            for reg in altpol_reg:(n_regs - 1)
                reg_dates[reg + 1] = reg_dates[reg] # shift all dates up by one
            end
            reg_dates[n_regs + 1] = reg_dates[n_regs] # add extra regime for new date
            reg_dates[altpol_reg] = altpol_date
        else
            altpol_reg = counter + 1 # just need to set the regime number
        end
    else
        altpol_reg = 2
        m <= Setting(:regime_dates, Dict{Int64, Date}(1 => date_presample_start(m), altpol_reg => altpol_date))
    end

    # Call setup_regime_switching_inds! to update indices
    setup_regime_switching_inds!(m; cond_type = cond_type)

    # Update regime_eqcond_info
    m <= Setting(:replace_eqcond, true)
    if haskey(get_settings(m), :regime_eqcond_info)
        reg_eqcond_info = get_setting(m, :regime_eqcond_info)
        reg_eqcond_info[altpol_reg] = EqcondEntry(altpol)
        if remove_extra_eqcond_entry # Remove any extra EqcondEntry values
            for reg in (altpol_reg + 1):get_setting(m, :n_regimes) # go to n_regimes b/c called setup_regime_switching_inds! already
                delete!(reg_eqcond_info, reg)
            end
        end
    else
        m <= Setting(:regime_eqcond_info, Dict{Int, EqcondEntry}(altpol_reg => EqcondEntry(altpol)))
    end

    return m
end
