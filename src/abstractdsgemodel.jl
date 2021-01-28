"""
```
AbstractDSGEModel{T} <: AbstractModel{T}
```

The AbstractDSGEModel is defined as a subtype of AbstractModel
"""
abstract type AbstractDSGEModel{T} <: ModelConstructors.AbstractModel{T} end

function Base.show(io::IO, m::AbstractDSGEModel)
    @printf io "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "no. states:             %i\n" n_states(m)
    @printf io "no. anticipated policy shocks: %i\n" n_mon_anticipated_shocks(m)
    @printf io "data vintage:           %s\n" data_vintage(m)
    @printf io "description:\n %s\n"          description(m)
end


# Retrieve model type
model_type(m::AbstractDSGEModel{T}) where T = T

# Number of anticipated policy shocks
n_anticipated_shocks(m::AbstractDSGEModel) = get_setting(m, :n_mon_anticipated_shocks)
n_anticipated_shocks_padding(m::AbstractDSGEModel) = get_setting(m, :n_mon_anticipated_shocks_padding)
n_mon_anticipated_shocks(m::AbstractDSGEModel) = get_setting(m, :n_mon_anticipated_shocks)
n_mon_anticipated_shocks_padding(m::AbstractDSGEModel) = get_setting(m, :n_mon_anticipated_shocks_padding)
n_anticipated_shocks_padding(m::AbstractDSGEModel) = get_setting(m, :n_mon_anticipated_shocks_padding)
n_z_anticipated_shocks(m::AbstractDSGEModel) = get_setting(m, :n_z_anticipated_shocks)
n_z_anticipated_shocks_padding(m::AbstractDSGEModel) = get_setting(m, :n_z_anticipated_shocks_padding)


# Dates, indices, number of periods for each regime
date_presample_start(m::AbstractDSGEModel) = get_setting(m, :date_presample_start)
date_mainsample_start(m::AbstractDSGEModel) = get_setting(m, :date_mainsample_start)
date_zlb_start(m::AbstractDSGEModel) = get_setting(m, :date_zlb_start)
date_zlb_end(m::AbstractDSGEModel) = get_setting(m, :date_zlb_end)
date_presample_end(m::AbstractDSGEModel) = Dates.lastdayofquarter(get_setting(m, :date_mainsample_start) - Dates.Month(3))
date_prezlb_end(m::AbstractDSGEModel) = Dates.lastdayofquarter(get_setting(m, :date_zlb_start) - Dates.Month(3))
date_mainsample_end(m::AbstractDSGEModel) = Dates.lastdayofquarter(get_setting(m, :date_forecast_start) - Dates.Month(3))
date_conditional_end(m::AbstractDSGEModel) = get_setting(m, :date_conditional_end)

index_presample_start(m::AbstractDSGEModel) = 1
index_mainsample_start(m::AbstractDSGEModel) = subtract_quarters(date_mainsample_start(m), date_presample_start(m)) + 1
index_zlb_start(m::AbstractDSGEModel) = subtract_quarters(date_zlb_start(m), date_presample_start(m)) + 1
index_forecast_start(m::AbstractDSGEModel) = subtract_quarters(date_forecast_start(m), date_presample_start(m)) + 1

"""
```
index_shockdec_start(m::AbstractDSGEModel)
```

Returns the index starting from which the shock decomposition is saved, where 1 is the index corresponding to date_mainsample_start(m).
"""
index_shockdec_start(m::AbstractDSGEModel) = subtract_quarters(date_shockdec_start(m), date_mainsample_start(m)) + 1

"""
```
index_shockdec_end(m::AbstractDSGEModel)
```

Returns the last index for which the shock decomposition is saved, where 1 is the index corresponding to date_mainsample_start(m).
"""
index_shockdec_end(m::AbstractDSGEModel) = subtract_quarters(date_shockdec_end(m), date_mainsample_start(m)) + 1

n_presample_periods(m::AbstractDSGEModel)   = subtract_quarters(date_mainsample_start(m), date_presample_start(m))
n_prezlb_periods(m::AbstractDSGEModel)      = subtract_quarters(date_zlb_start(m), date_mainsample_start(m))
n_zlb_periods(m::AbstractDSGEModel)         = subtract_quarters(date_forecast_start(m), date_zlb_start(m))
n_mainsample_periods(m::AbstractDSGEModel)  = subtract_quarters(date_forecast_start(m), date_mainsample_start(m))
n_conditional_periods(m::AbstractDSGEModel) = subtract_quarters(date_conditional_end(m), date_mainsample_end(m))

inds_presample_periods(m::AbstractDSGEModel)  = collect(index_presample_start(m):(index_mainsample_start(m)-1))
inds_prezlb_periods(m::AbstractDSGEModel)     = collect(index_mainsample_start(m):(index_zlb_start(m)-1))
inds_zlb_periods(m::AbstractDSGEModel)        = collect(index_zlb_start(m):(index_forecast_start(m)-1))
inds_mainsample_periods(m::AbstractDSGEModel) = collect(index_mainsample_start(m):(index_forecast_start(m)-1))

# Convenience functions for working with heterogeneous agent models
# that differentiate between backward looking "state" variables and "jump" variables
n_backward_looking_states(m::AbstractDSGEModel) = get_setting(m, :n_backward_looking_states)
n_jumps(m::AbstractDSGEModel) = get_setting(m, :n_jumps)
n_model_states(m::AbstractDSGEModel) = get_setting(m, :n_model_states)

# The numbers for n_states, and n_jumps assumes normalization
# There is a procedure in klein_solve that normalizes the state variable grids
# by removing an entry from them. Hence the number of states and jumps being
# tracked should be 1 less if normalized
# However, the number of jumps and states unnormalized are required for
# the construction of the Jacobian, hence the reason for these helpers
n_jumps_unnormalized(m::AbstractDSGEModel) = n_jumps(m) + get_setting(m, :jumps_normalization_factor)
n_backward_looking_states_unnormalized(m::AbstractDSGEModel) = n_backward_looking_states(m) + get_setting(m, :backward_looking_states_normalization_factor)
n_model_states_unnormalized(m::AbstractDSGEModel) = n_jumps_unnormalized(m) + n_backward_looking_states_unnormalized(m)
n_model_states_original(m::AbstractDSGEModel) = get_setting(m, :n_model_states_original)

"""
```
AbstractCTModel{T} <: AbstractDSGEModel{T}
```

The AbstractCTModel is defined as a subtype of AbstractDSGEModel to accomodate the
numerical methods and procedures specific to continuous time models.
"""
abstract type AbstractCTModel{T} <: AbstractDSGEModel{T} end

n_states(m::AbstractCTModel) = sum(map(i -> length(collect(m.endogenous_states)[i][2]), 1:length(keys(m.endogenous_states))))
n_shocks_expectational(m::AbstractCTModel) = sum(map(i -> length(collect(m.expected_shocks)[i][2]), 1:length(keys(m.expected_shocks))))

"""
```
AbstractHetModel{T} <: AbstractDSGEModel{T}
```

The AbstractHetModel is defined as a subtype of AbstractDSGEModel to accomodate a bunch of stuff, but for now just different impulse response functions.
"""
abstract type AbstractHetModel{T} <: AbstractDSGEModel{T} end

"""
```
AbstractRepModel{T} <: AbstractDSGEModel{T}
```

The AbstractRepresentativeModel is defined as a subtype of AbstractDSGEModel to accomodate a bunch of stuff, but for now just different impulse response functions.
"""
abstract type AbstractRepModel{T} <: AbstractDSGEModel{T} end

"""
```
get_dict(m, class, index)
```
"""
function get_dict(m::AbstractDSGEModel, class::Symbol)
    if class == :states
        m.endogenous_states
    elseif class == :obs
        m.observables
    elseif class == :pseudo
        m.pseudo_observables
    elseif class in [:shocks, :stdshocks]
        m.exogenous_shocks
    else
        throw(ArgumentError("Invalid class: $class. Must be :states, :obs, :pseudo, :shocks, or :stdshocks"))
    end
end

"""
```
get_key(m, class, index)
```
Returns the name of the state (`class = :states`), observable (`:obs`),
pseudo-observable (`:pseudo`), or shock (`:shocks` or `:stdshocks`)
corresponding to the given `index`.
"""
function get_key(m::AbstractDSGEModel, class::Symbol, index::Int)
    dict = get_dict(m, class)
    out = Base.filter(key -> dict[key] == index, collect(keys(dict)))
    if length(out) == 0
        error("Key corresponding to index $index not found for class: $class")
    elseif length(out) > 1
        error("Multiple keys corresponding to index $index found for class: $class")
    else
        return out[1]
    end
end


# Parse population mnemonic into 2 Nullable{Symbol}s from one
function parse_population_mnemonic(m::AbstractDSGEModel)
    mnemonic = get_setting(m, :population_mnemonic)
    if isnull(mnemonic)
        return [Nullable{Symbol}(), Nullable{Symbol}()]
    else
        return map(s -> Nullable(Symbol(s)), split(string(get(mnemonic)), DSGE_DATASERIES_DELIM))
    end
end

# From an augmented state space with anticipated policy shocks, get indices
# corresponding to pre-ZLB states, shocks, and observables
function inds_states_no_ant(m::AbstractDSGEModel)
    if n_mon_anticipated_shocks(m) > 0
        ind_ant1 = m.endogenous_states[:rm_tl1]
        ind_antn = m.endogenous_states[Symbol("rm_tl$(n_mon_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_states_augmented(m)]
    else
        return collect(1:n_states_augmented(m))
    end
end

function inds_shocks_no_ant(m::AbstractDSGEModel)
    if n_mon_anticipated_shocks(m) > 0
        ind_ant1 = m.exogenous_shocks[:rm_shl1]
        ind_antn = m.exogenous_shocks[Symbol("rm_shl$(n_mon_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_shocks_exogenous(m)]
    else
        return collect(1:n_shocks_exogenous(m))
    end
end

function inds_obs_no_ant(m::AbstractDSGEModel)
    if n_mon_anticipated_shocks(m) > 0
        ind_ant1 = m.observables[:obs_nominalrate1]
        ind_antn = m.observables[Symbol("obs_nominalrate$(n_mon_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_observables(m)]
    else
        return collect(1:n_observables(m))
    end
end

# From an augmented state space with integrated series (e.g. unit root),
# get indices corresponding to stationary states
function inds_states_no_integ_series(m::AbstractDSGEModel)
    if haskey(get_settings(m), :integrated_series)
        inds = map(i -> m.endogenous_states_augmented[i],
                   get_setting(m, :integrated_series))
        return setdiff(1:n_states_augmented(m), inds)
    else
        return collect(1:n_states_augmented(m))
    end
end

# Interface for accessing parameters
get_parameters(m::AbstractDSGEModel) = m.parameters

# Interface for accessing rng
get_rng(m::AbstractDSGEModel) = m.rng

# Interface for accessing settings dictionary
get_settings(m::AbstractDSGEModel) = hasproperty(m, :testing) ? (m.testing ? m.test_settings : m.settings) : m.settings

# Interface for accessing observables dictionary
get_observables(m::AbstractDSGEModel) = m.observables
get_pseudo_observables(m::AbstractDSGEModel) = m.pseudo_observables
get_exogenous_shocks(m::AbstractDSGEModel) = m.exogenous_shocks

# Interface for data
cond_vintage(m::AbstractDSGEModel)    = get_setting(m, :cond_vintage)
cond_id(m::AbstractDSGEModel)         = get_setting(m, :cond_id)
cond_full_names(m::AbstractDSGEModel) = get_setting(m, :cond_full_names)
cond_semi_names(m::AbstractDSGEModel) = get_setting(m, :cond_semi_names)
use_population_forecast(m::AbstractDSGEModel) = get_setting(m, :use_population_forecast)
hpfilter_population(m::AbstractDSGEModel)     = get_setting(m, :hpfilter_population)

# Interface for general computation settings
use_parallel_workers(m::AbstractDSGEModel)    = get_setting(m, :use_parallel_workers)

# Interface for estimation settings
reoptimize(m::AbstractDSGEModel)          = get_setting(m, :reoptimize)
calculate_hessian(m::AbstractDSGEModel) = get_setting(m, :calculate_hessian)
hessian_path(m::AbstractDSGEModel)      = get_setting(m, :hessian_path)
n_hessian_test_params(m::AbstractDSGEModel) = get_setting(m, :n_hessian_test_params)

# Interface for Metropolis-Hastings settings
n_mh_blocks(m::AbstractDSGEModel)       =  get_setting(m, :n_mh_blocks)
n_mh_param_blocks(m::AbstractDSGEModel) =  get_setting(m, :n_mh_param_blocks)
n_mh_simulations(m::AbstractDSGEModel)  =  get_setting(m, :n_mh_simulations)
n_mh_burn(m::AbstractDSGEModel)         =  get_setting(m, :n_mh_burn)
mh_thin(m::AbstractDSGEModel)           =  get_setting(m, :mh_thin)

# Interface for forecast settings
date_forecast_start(m::AbstractDSGEModel)   = get_setting(m, :date_forecast_start)
forecast_block_size(m::AbstractDSGEModel)   = get_setting(m, :forecast_block_size)
forecast_start_block(m::AbstractDSGEModel)  = get_setting(m, :forecast_start_block)
forecast_input_file_overrides(m::AbstractDSGEModel) = get_setting(m, :forecast_input_file_overrides)
forecast_uncertainty_override(m::AbstractDSGEModel) = get_setting(m, :forecast_uncertainty_override)
forecast_smoother(m::AbstractDSGEModel)     = get_setting(m, :forecast_smoother)
forecast_tdist_df_val(m::AbstractDSGEModel) = get_setting(m, :forecast_tdist_df_val)
forecast_tdist_shocks(m::AbstractDSGEModel) = get_setting(m, :forecast_tdist_shocks)
forecast_zlb_value(m::AbstractDSGEModel)    = get_setting(m, :forecast_zlb_value)
impulse_response_horizons(m::AbstractDSGEModel) = get_setting(m, :impulse_response_horizons)
n_shockdec_periods(m::AbstractDSGEModel)    = index_shockdec_end(m) - index_shockdec_start(m) + 1

# Function for getting the alternative policy. Note that, in the case of imperfect awareness,
# i.e. when uncertain_altpolicy is true, this function will return the policy that is actually implemented,
# instead of other potential alternative policies which households believe could occur.
#
# Furthermore, do NOT delete the backup cases when regime_eqcond_info is not defined, namely
# the check for whether :alternative_policy is a Setting. This last check is necessary for
# the scenarios code to continue working.
alternative_policy(m::AbstractDSGEModel) = haskey(get_settings(m), :regime_eqcond_info) && haskey(get_settings(m), :n_regimes) &&
    haskey(get_setting(m, :regime_eqcond_info), get_setting(m, :n_regimes)) ?
    get_setting(m, :regime_eqcond_info)[get_setting(m, :n_regimes)].alternative_policy :
    (haskey(get_settings(m), :alternative_policy) ? get_setting(m, :alternative_policy) : AltPolicy(:historical, eqcond, solve))

# Some additional date settings related to forecasts
function date_forecast_end(m::AbstractDSGEModel)
    if haskey(get_settings(m), :date_forecast_end)
        return get_setting(m, :date_forecast_end)
    else
        date = date_forecast_start(m) + Dates.Month(3 * (forecast_horizons(m)-1))
        return Dates.lastdayofquarter(date)
    end
end

function forecast_horizons(m::AbstractDSGEModel; cond_type::Symbol = :none)
    horizons = get_setting(m, :forecast_horizons)
    if cond_type == :none
        return horizons
    else
        return horizons - n_conditional_periods(m)
    end
end

function date_shockdec_start(m::AbstractDSGEModel)
    startdate = get_setting(m, :shockdec_startdate)
    if (Nullables.isnull(startdate)) | (startdate===nothing)
        return DSGE.date_mainsample_start(m)
    else
        return get(startdate)
    end
end

function date_shockdec_end(m::AbstractDSGEModel)
    enddate =  get_setting(m, :shockdec_enddate)
    if !isnull(enddate)
        return get(enddate)
    else
        return date_forecast_end(m)
    end
end

"""
```
load_parameters_from_file(m::AbstractDSGEModel,path::String)
```
Returns a vector of parameters, read from a file, suitable for updating `m`.
"""
function load_parameters_from_file(m::AbstractDSGEModel, path::String)

    if isfile(path) && splitext(path)[2] == ".h5"
        x  = h5open(path, "r") do file
            try
                read(file, "params")
            catch
                error("$path does not contain variable params")
            end
        end
    else
        error("$path is not a valid HDF5 file.")
    end

    if m.spec=="smets_wouters" && length(x)==42
        x = vcat(x, zeros(7))
    end

    @assert length(x) == length(m.parameters)
    @assert eltype(x) == typeof(m.parameters[1].value)
    return x
end

"""
```
specify_mode!(m::AbstractDSGEModel, mode_file::String=""; verbose=:low)
```

Updates the values of `m.parameters` with the values from
`mode_file`. Sets `reoptimize` setting to `false`.

Usage: should be run before calling `estimate(m)`, e.g.:

    m = Model990()
    specify_mode!(m, modefile)
    estimate(m)
"""
function specify_mode!(m::AbstractDSGEModel, mode_file::String = ""; verbose=:low)

    m <= Setting(:reoptimize, false)

    if mode_file == ""
        mode_file = inpath(m, "user", "paramsmode.h5")
    end
    mode_file = normpath(mode_file)

    DSGE.update!(m,load_parameters_from_file(m,mode_file))

    println(verbose, :low, "Loaded previous mode from $mode_file.")
end

"""
```
specify_hessian!(m::AbstractDSGEModel, path::String=""; verbose=:low)
```

Specify a Hessian matrix calculated at the posterior mode to use in the model estimation. If
no path is provided, will attempt to detect location.
"""
function specify_hessian!(m::AbstractDSGEModel, path::String=""; verbose=:low)
    if isempty(path)
        path = inpath(m, "user", "hessian.h5")
    end

    if isfile(path) && splitext(path)[2] == ".h5"
        m <= Setting(:hessian_path, normpath(abspath(path)))
    else
        error("Invalid input Hessian file: $path",)
    end

    m <= Setting(:calculate_hessian, false)

    println(verbose, :low, "Specified hessian from $path.")
end


"""
```
transform_to_model_space!(m::AbstractDSGEModel, values::Vector{T}; regime_switching::Bool = false) where T<:AbstractFloat
```

Transforms `values` from the real line to the model space, and assigns `values[i]` to
`m.parameters[i].value` for non-steady-state parameters. Recomputes the steady-state
parameter values.

### Arguments
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.
- `regime_switching`: set to true if the model's parameters are regime-switching
"""
function transform_to_model_space!(m::AbstractDSGEModel, values::Vector{T};
                                   regime_switching::Bool = false) where {T<:AbstractFloat}
    new_values = transform_to_model_space(m.parameters, values; regime_switching = regime_switching)
    DSGE.update!(m, new_values)
    steadystate!(m)
end

"""
```
update!(m::AbstractDSGEModel, values::Vector{T}) where T<:AbstractFloat
```

Update `m.parameters` with `values`, recomputing the steady-state parameter values.

### Arguments:
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.
"""
function update!(m::AbstractDSGEModel, values::Vector{T}) where {T <: Real}
    ModelConstructors.update!(m.parameters, values)
    steadystate!(m)
end

"""
```
update!(m::AbstractDSGEModel, values::ParameterVector{T};
    regime_switching::Bool = false, toggle::Bool = true) where T
```
Update `m.parameters` with `values`, recomputing the steady-state parameter values.

### Arguments
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.

### Keyword
- `regime_switching`: if true, then we assume the parameters are regime-switching,
    in which case `update!` assumes the `value` field
    of each parameter in values` holds the parameter value in the first regime, and
    then we update the field `regimes` for each parameter
- `toggle`: if true, we call `ModelConstructors.toggle_regime!(values)` before
    updating any values to ensure the `value` field of the parameters in `values`
    correspond to regime 1 values.
"""
function update!(m::AbstractDSGEModel, values::ParameterVector{T};
                 regime_switching::Bool = false, toggle::Bool = true) where {T <: Real}

    # Update regime-switching if length of `values` exceeds m.parameters
    if regime_switching
        if toggle
            ModelConstructors.toggle_regime!(values, 1)
        end

        # Update first-regime values
        ModelConstructors.update!(m.parameters, [θ.value for θ in values])

        # Update remaining regimes
        for (i, para) in enumerate(m.parameters)
            if !isempty(para.regimes)
                for (ind, val) in para.regimes[:value]
                    if ind == 1
                        ModelConstructors.set_regime_val!(para, 1, para.value)
                    else
                        ModelConstructors.set_regime_val!(para, ind, regime_val(values[i], ind))
                    end
                end
            end
        end
    else
        ModelConstructors.update!(m.parameters, [θ.value for θ in values])
    end

    steadystate!(m)
end

"""
```
mutable struct ShockGroup
```

The `ShockGroup` mutable struct is used in `prepare_means_table_shockdec` and
`plot_shock_decompositions`. When plotting shock decompositions, we usually want
to group the shocks into categories (financial, monetary policy, etc.) so that
the resulting grouped bar plot is legible.

### Fields

- `name::String`
- `shocks::Vector{Symbol}`
- `color::Colorant`
"""
mutable struct ShockGroup
    name::String
    shocks::Vector{Symbol}
    color::Colorant
end

function ShockGroup(name::String, shocks::Vector{Symbol}, color_name::Symbol)
    color = parse(Colorant, color_name)
    return ShockGroup(name, shocks, color)
end

mutable struct SteadyStateConvergenceError <: Exception
    msg::String
end
SteadyStateConvergenceError() = SteadyStateConvergenceError("SteadyState didn't converge")
Base.showerror(io::IO, ex::SteadyStateConvergenceError) = print(io, ex.msg)


"""
```
setup_regime_switching_inds!(m::AbstractDSGEModel; cond_type::Symbol = :none,
    temp_altpolicy_in_cond_regimes::Bool = false)
```

calculates the indices needed to solve and forecast a model with regime-switching.

### Keywords
- `cond_type`: the correct regime indices for forecasting depend on whether the forecast is conditional or not
- `temp_altpolicy_in_cond_regimes`: If true, then temporary alternative policies can occur
    during the conditional forecast horizon, and regime indices (namely `n_rule_periods`) will be calculated accordingly.
"""
function setup_regime_switching_inds!(m::AbstractDSGEModel; cond_type::Symbol = :none,
                                      temp_altpolicy_in_cond_regimes::Bool = false)

    if haskey(get_settings(m), :regime_switching) ? !get_setting(m, :regime_switching) : true
        @warn "The setting :regime_switching is either false or is not defined yet. Updating the setting to be true."
        m <= Setting(:regime_switching, true)
    end

    n_hist_regimes = 0
    n_cond_regimes = 0
    n_regimes = length(get_setting(m, :regime_dates))
    post_cond_end = iterate_quarters(date_conditional_end(m), 1) # Period after conditional forecasting ends
    set_reg_forecast_start   = false # These flags are needed to tell whether or not we actually set these regimes.
    set_post_conditional_end = false # W/out these flags, if m already has reg_forecast_start, then it won't be properly set
    reg_forecast_start_str = "Regime in which the forecast starts."
    reg_post_conditional_end_str = "Regime one period after the conditional forecast ends."
    for (key, val) in get_setting(m, :regime_dates)
        if val == date_forecast_start(m)
            m <= Setting(:reg_forecast_start, key, reg_forecast_start_str)
            set_reg_forecast_start = true
        end
        if val == post_cond_end
            m <= Setting(:reg_post_conditional_end, key, reg_post_conditional_end_str)
            set_post_conditional_end = true
        end
        if val < date_forecast_start(m)
            n_hist_regimes += 1
        end
        if date_forecast_start(m) <= val <= date_conditional_end(m)
            n_cond_regimes += 1
        end
    end
    if !set_reg_forecast_start
        # Then the forecast begins in the middle of a regime
        m <= Setting(:reg_forecast_start, findlast(sort!(collect(values(get_setting(m, :regime_dates)))) .< date_forecast_start(m)),
                     reg_forecast_start_str)
    end
    if !set_post_conditional_end
        # Then the conditional forecast ends in the middle of a regime
        m <= Setting(:reg_post_conditional_end,
                     findlast(sort!(collect(values(get_setting(m, :regime_dates)))) .<= date_conditional_end(m)),
                     reg_post_conditional_end_str) # or .< post_cond_end
    end

    # Infer number of regimes in the forecast horizon from reg_forecast_start
    n_fcast_regimes = n_regimes - get_setting(m, :reg_forecast_start) + 1

    m <= Setting(:n_regimes, n_regimes, "Total number of regimes")
    m <= Setting(:n_hist_regimes, n_hist_regimes, "Number of regimes in the history")
    m <= Setting(:n_fcast_regimes, n_fcast_regimes, "Number of regimes in the forecast horizon")
    m <= Setting(:n_cond_regimes, n_cond_regimes, "Number of regime switches during the conditional forecast horizon")

    # Number of periods that the rule is in place is (normally)
    # n_regimes - n_hist_regimes - 1 (-1 b/c in last regime, go back to normal rule)
    # If conditional forecasting occurs and the rule occurs afterward, we need to subtract n_cond_regimes.
    # If rule occurs during the conditional forecast horizon or before, then we need to add
    # extra periods to n_rule_periods
    # TODO: Update this temp_altpolicy_in_cond_regimes to reflect regime_eqcond_info (may no longer be needed)
    cond_adj = if temp_altpolicy_in_cond_regimes
        haskey(get_settings(m), :gensys2_first_regime) ?
            (get_setting(m, :gensys2_first_regime) - get_setting(m, :n_hist_regimes) - 1) : 0
    else
        get_setting(m, :n_cond_regimes)
    end
    m <= Setting(:n_rule_periods, (cond_type == :none) ? n_regimes - (get_setting(m, :n_hist_regimes) + 1) :
                 n_regimes - (get_setting(m, :n_hist_regimes) + 1 + cond_adj),
                 "Number of periods during which the (temporary) alternative policy applies")
    return m
end

# Dummy function for eqcond, relevant when using gensys2
function eqcond(m::AbstractDSGEModel, Γ0::AbstractMatrix{S}, Γ1::AbstractMatrix{S},
                C::AbstractVector{S}, Ψ::AbstractMatrix{S}, Π::AbstractMatrix{S}) where {S <: Real}
    return Γ0, Γ1, C, Ψ, Π
end

"""
```
setup_param_regimes!(m::AbstractDSGEModel, param_regs::Matrix{Int} = []

Function to set up the model with parameter regime switching
for estimation. param_regs should be a matrix where for each row r,
the i{th} element is the parameter regime for the i{th} model regime
for parameter r (including parameters with no regime-switching: set
the rows for non-regime-switching parameters to all 1s, although any
value will give the same result).
```
"""
function setup_param_regimes!(m::AbstractDSGEModel, param_mat::Array{Int, 2} = Matrix{Int}(undef, 0, 0))

    param_reg   = Dict{Symbol, Dict{Int, Int}}()
    nmodel_regs = haskey(get_settings(m), :n_regimes) ? get_setting(m, :n_regimes) : 1

    if isempty(param_mat)
        param_mat = ones(n_parameters(m), nmodel_regs)
    end

    @assert nmodel_regs == size(param_mat, 2) "The number of columns in `param_mat` must match the number of model regimes"

    for i in 1:n_parameters(m)
        reg_dict = Dict{Int, Int}()
        for reg in 1:nmodel_regs
            reg_dict[reg] = param_mat[i, reg]
        end
        param_reg[m.parameters[i].key] = reg_dict
    end

    m <= Setting(:model2para_regime, param_reg)

    return m
end
