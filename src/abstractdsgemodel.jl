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
    @printf io "no. anticipated shocks: %i\n" n_anticipated_shocks(m)
    @printf io "data vintage:           %s\n" data_vintage(m)
    @printf io "description:\n %s\n"          description(m)
end

# Number of anticipated policy shocks
n_anticipated_shocks(m::AbstractDSGEModel) = get_setting(m, :n_anticipated_shocks)
n_anticipated_shocks_padding(m::AbstractDSGEModel) = get_setting(m, :n_anticipated_shocks_padding)

# Dates, indices, number of periods for each regime
date_presample_start(m::AbstractDSGEModel) = get_setting(m, :date_presample_start)
date_mainsample_start(m::AbstractDSGEModel) = get_setting(m, :date_mainsample_start)
date_zlb_start(m::AbstractDSGEModel) = get_setting(m, :date_zlb_start)

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
AbstractRepModel{T} <: AbstractDSGEModel{T}
```

The AbstractRepresentativeModel is defined as a subtype of AbstractDSGEModel to accomodate a bunch of stuff, but for now just different impulse response functions.
"""
abstract type AbstractRepModel{T} <: AbstractDSGEModel{T} end


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
    if n_anticipated_shocks(m) > 0
        ind_ant1 = m.endogenous_states[:rm_tl1]
        ind_antn = m.endogenous_states[Symbol("rm_tl$(n_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_states_augmented(m)]
    else
        return collect(1:n_states_augmented(m))
    end
end

function inds_shocks_no_ant(m::AbstractDSGEModel)
    if n_anticipated_shocks(m) > 0
        ind_ant1 = m.exogenous_shocks[:rm_shl1]
        ind_antn = m.exogenous_shocks[Symbol("rm_shl$(n_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_shocks_exogenous(m)]
    else
        return collect(1:n_shocks_exogenous(m))
    end
end

function inds_obs_no_ant(m::AbstractDSGEModel)
    if n_anticipated_shocks(m) > 0
        ind_ant1 = m.observables[:obs_nominalrate1]
        ind_antn = m.observables[Symbol("obs_nominalrate$(n_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_observables(m)]
    else
        return collect(1:n_observables(m))
    end
end

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

# Interface for alternative policy settings
alternative_policy(m::AbstractDSGEModel) = get_setting(m, :alternative_policy)

function date_forecast_end(m::AbstractDSGEModel)
    date = date_forecast_start(m) + Dates.Month(3 * (forecast_horizons(m)-1))
    return Dates.lastdayofquarter(date)
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
transform_to_model_space!(m::AbstractDSGEModel, values::Vector{T}) where T<:AbstractFloat
```

Transforms `values` from the real line to the model space, and assigns `values[i]` to
`m.parameters[i].value` for non-steady-state parameters. Recomputes the steady-state
paramter values.

### Arguments
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.
"""
function transform_to_model_space!(m::AbstractDSGEModel, values::Vector{T}) where {T<:AbstractFloat}
    new_values = transform_to_model_space(m.parameters, values)
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
function update!(m::AbstractDSGEModel, values::Vector{T}) where T<:AbstractFloat
    ModelConstructors.update!(m.parameters, values)
    steadystate!(m)
end

"""
```
update!(m::AbstractDSGEModel, values::ParameterVector{T}) where T
```
Update `m.parameters` with `values`, recomputing the steady-state parameter values.
### Arguments:
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.
"""
function update!(m::AbstractDSGEModel, values::ParameterVector{T}) where T
    ModelConstructors.update!(m.parameters, [θ.value for θ in values])
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
