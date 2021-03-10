abstract type AbstractAltPolicy end

Base.string(a::AbstractAltPolicy) = string(a.key)

"""
```
mutable struct AltPolicy
```

Types defining an alternative policy rule.

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
mutable struct AltPolicy <: AbstractAltPolicy
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

"""
```
get_altpolicy(key::Symbol)
```

returns the `AltPolicy` for a alternative policy identified by `key`.
"""
function get_altpolicy(key::Symbol)
    altpol = if key == :flexible_ait
        flexible_ait()
    elseif key == :zero_rate
        zero_rate()
    elseif key == :zlb_rule
        zlb_rule()
    elseif key == :default_policy
        default_policy()
    elseif key == :taylor_rule
        taylor_rule()
    elseif key == :smooth_ait_gdp_alt
        smooth_ait_gdp_alt()
    elseif key == :smooth_ait_gdp
        smooth_ait_gdp()
    elseif key == :ait
        ait()
    elseif key == :ngdp
        ngdp()
    elseif key == :rw_zero_rate
        rw_zero_rate()
    elseif key == :rw
        rw()
    else
        error("AltPolicy $(key) has not been added to `get_altpolicy` as having a method for replacing equilibrium conditions")
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
mutable struct MultiPeriodAltPolicy
```

Types defining an alternative policy rule.

### Fields

- `key::Symbol`: alternative policy identifier

- `regime_eqcond_info::AbstractDict{Int64, EqcondEntry}`: a dictionary mapping
    regimes to equilibrium conditions which replace the default ones in a
    given regime.

- `gensys2::Bool`: if true, the multi-period alternative policy needs
    to call `gensys2` instead of `gensys` to work.

- `infoset::Union{Vector{UnitRange{Int64}}, Nothing}`: either a vector specifying
    the information set used for expectations in the measurement equation or
    `nothing` to indicate myopia in expectations across regimes.

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
mutable struct MultiPeriodAltPolicy{S <: AbstractDict{Int64, EqcondEntry}} <: AbstractAltPolicy
    key::Symbol
    n_regimes::Int64
    regime_eqcond_info::S
    gensys2::Bool
    temporary_altpolicy_names::Union{Vector{Symbol}, Nothing}
    temporary_altpolicy_length::Int64
    infoset::Union{Vector{UnitRange{Int64}}, Nothing}
    perfect_cred_regime_mapping::Union{Dict{Int64, Int64}, Nothing}
    forecast_init::Function
    color::Colorant
    linestyle::Symbol
end

function MultiPeriodAltPolicy(key::Symbol, n_regimes::Int64, regime_eqcond_info::AbstractDict{Int64, EqcondEntry};
                              gensys2::Bool = false,
                              temporary_altpolicy_names::Union{Vector{Symbol}, Nothing} = nothing,
                              temporary_altpolicy_length::Int =  # default assumes the last regime is a
                              length(regime_eqcond_info) - 1, # lift-off regime, and all others are temporary ones.
                              infoset::Union{Vector{UnitRange{Int64}}, Nothing} = nothing,
                              perfect_cred_regime_mapping::Union{Dict{Int64, Int64}, Nothing} = nothing,
                              forecast_init::Function = identity,
                              color::Colorant = RGB(0., 0., 1.),
                              linestyle::Symbol = :solid)

    MultiPeriodAltPolicy(key, n_regimes, regime_eqcond_info, gensys2, temporary_altpolicy_names,
                         temporary_altpolicy_length, infoset,
                         perfect_cred_regime_mapping, forecast_init, color, linestyle)
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
