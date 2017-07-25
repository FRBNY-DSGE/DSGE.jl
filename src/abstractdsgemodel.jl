abstract AbstractModel{T}

function Base.show(io::IO, m::AbstractModel)
    @printf io "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "no. states:             %i\n" n_states(m)
    @printf io "no. anticipated shocks: %i\n" n_anticipated_shocks(m)
    @printf io "data vintage:           %s\n" data_vintage(m)
    @printf io "description:\n %s\n"          description(m)
end

@inline function Base.getindex(m::AbstractModel, i::Integer)
    if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]
    end
end

# need to define like this so we can disable bounds checking
@inline function Base.getindex(m::AbstractModel, k::Symbol)
    i = m.keys[k]
    @inbounds if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]
    end
end

@inline function Base.setindex!(m::AbstractModel, value::Number, i::Integer)
    if i <= (j = length(m.parameters))
        param = m.parameters[i]
        param.value = value
        if isa(param, ScaledParameter)
            param.scaledvalue = param.scaling(value)
        end
        return param
    else
        steady_state_param = m.steady_state[i-j]
        steady_state_param.value = value
        return steady_state_param
    end
end

"""
```
setindex!(m::AbstractModel, param::AbstractParameter, i::Integer)
```

If `i`<length(m.parameters), overwrites m.parameters[i] with
param. Otherwise, overwrites m.steady_state[i-length(m.parameters).
"""
@inline function Base.setindex!(m::AbstractModel, param::AbstractParameter, i::Integer)
    if i <= (j = length(m.parameters))
        m.parameters[i] = param
    else
        m.steady_state[i-j] = param
    end
    return param
end

Base.setindex!(m::AbstractModel, value, k::Symbol) = Base.setindex!(m, value, m.keys[k])


"""
```
(<=){T}(m::AbstractModel{T}, p::AbstractParameter{T})
```

Syntax for adding a parameter to a model: m <= parameter.
NOTE: If `p` is added to `m` and length(m.steady_state) > 0, `keys(m)` will not generate the
index of `p` in `m.parameters`.
"""
function (<=){T}(m::AbstractModel{T}, p::AbstractParameter{T})

    if !in(p.key, keys(m.keys))

        new_param_index = length(m.keys) + 1

        # grow parameters and add the parameter
        push!(m.parameters, p)

        # add parameter location to dict
        setindex!(m.keys, new_param_index, p.key)
    else
        # overwrite the previous parameter with the new one
        setindex!(m, p, p.key)
    end
end


"""
```
(<=){T}(m::AbstractModel{T}, ssp::SteadyStateParameter)
```

Add a new steady-state value to the model by appending `ssp` to the `m.steady_state` and
adding `ssp.key` to `m.keys`.
"""
function (<=){T}(m::AbstractModel{T}, ssp::SteadyStateParameter)

    if !in(ssp.key, keys(m.keys))
        new_param_index = length(m.keys) + 1

        # append ssp to steady_state vector
        push!(m.steady_state, ssp)

        # add parameter location to dict
        setindex!(m.keys, new_param_index, ssp.key)
    else
        # overwrite the previous parameter with the new one
        setindex!(m, ssp, ssp.key)
    end
end

Distributions.logpdf(m::AbstractModel) = logpdf(m.parameters)
Distributions.pdf(m::AbstractModel) = exp(logpdf(m))

# Number of anticipated policy shocks
n_anticipated_shocks(m::AbstractModel) = get_setting(m, :n_anticipated_shocks)
n_anticipated_shocks_padding(m::AbstractModel) = get_setting(m, :n_anticipated_shocks_padding)

# Dates, indices, number of periods for each regime
date_presample_start(m::AbstractModel) = get_setting(m, :date_presample_start)
date_mainsample_start(m::AbstractModel) = get_setting(m, :date_mainsample_start)
date_zlb_start(m::AbstractModel) = get_setting(m, :date_zlb_start)

date_presample_end(m::AbstractModel) = Dates.lastdayofquarter(get_setting(m, :date_mainsample_start) - Dates.Month(3))
date_prezlb_end(m::AbstractModel) = Dates.lastdayofquarter(get_setting(m, :date_zlb_start) - Dates.Month(3))
date_mainsample_end(m::AbstractModel) = Dates.lastdayofquarter(get_setting(m, :date_forecast_start) - Dates.Month(3))
date_conditional_end(m::AbstractModel) = get_setting(m, :date_conditional_end)

index_presample_start(m::AbstractModel) = 1
index_mainsample_start(m::AbstractModel) = subtract_quarters(date_mainsample_start(m), date_presample_start(m)) + 1
index_zlb_start(m::AbstractModel) = subtract_quarters(date_zlb_start(m), date_presample_start(m)) + 1
index_forecast_start(m::AbstractModel) = subtract_quarters(date_forecast_start(m), date_presample_start(m)) + 1

"""
```
index_shockdec_start(m::AbstractModel)
```

Returns the index starting from which the shock decomposition is saved, where 1 is the index corresponding to date_mainsample_start(m).
"""
index_shockdec_start(m::AbstractModel) = subtract_quarters(date_shockdec_start(m), date_mainsample_start(m)) + 1

"""
```
index_shockdec_end(m::AbstractModel)
```

Returns the last index for which the shock decomposition is saved, where 1 is the index corresponding to date_mainsample_start(m).
"""
index_shockdec_end(m::AbstractModel) = subtract_quarters(date_shockdec_end(m), date_mainsample_start(m)) + 1

n_presample_periods(m::AbstractModel)   = subtract_quarters(date_mainsample_start(m), date_presample_start(m))
n_prezlb_periods(m::AbstractModel)      = subtract_quarters(date_zlb_start(m), date_mainsample_start(m))
n_zlb_periods(m::AbstractModel)         = subtract_quarters(date_forecast_start(m), date_zlb_start(m))
n_mainsample_periods(m::AbstractModel)  = subtract_quarters(date_forecast_start(m), date_mainsample_start(m))
n_conditional_periods(m::AbstractModel) = subtract_quarters(date_conditional_end(m), date_mainsample_end(m))

inds_presample_periods(m::AbstractModel) = collect(index_presample_start(m):(index_mainsample_start(m)-1))
inds_prezlb_periods(m::AbstractModel) = collect(index_mainsample_start(m):(index_zlb_start(m)-1))
inds_zlb_periods(m::AbstractModel) = collect(index_zlb_start(m):(index_forecast_start(m)-1))
inds_mainsample_periods(m::AbstractModel) = collect(index_mainsample_start(m):(index_forecast_start(m)-1))

# Number of a few things that are useful
n_states(m::AbstractModel)                 = length(m.endogenous_states)
n_states_augmented(m::AbstractModel)       = n_states(m) + length(m.endogenous_states_augmented)
n_shocks_exogenous(m::AbstractModel)       = length(m.exogenous_shocks)
n_shocks_expectational(m::AbstractModel)   = length(m.expected_shocks)
n_equilibrium_conditions(m::AbstractModel) = length(m.equilibrium_conditions)
n_observables(m::AbstractModel)            = length(m.observables)
n_parameters(m::AbstractModel)             = length(m.parameters)
n_parameters_steady_state(m::AbstractModel)= length(m.steady_state)
n_parameters_free(m::AbstractModel)        = sum([!α.fixed for α in m.parameters])

function n_pseudoobservables(m::AbstractModel)
    if forecast_pseudoobservables(m)
        pseudo, _ = pseudo_measurement(m)
        return length(pseudo)
    else
        return 0
    end
end

"""
```
get_key(m, class, index)
```

Returns the name of the state (`class = :state`), observable (`obs`),
pseudo-observable (`pseudo`), or shock (`shock`) corresponding to the given
`index`.
"""
function get_key(m::AbstractModel, class::Symbol, index::Int)
    dict = if class == :state
        m.endogenous_states
    elseif class == :obs
        m.observables
    elseif class == :pseudo
        _, pseudo_mapping = pseudo_measurement(m)
        pseudo_mapping.inds
    elseif class in [:shock, :stdshock]
        m.exogenous_shocks
    else
        throw(ArgumentError("Invalid class :$class. Must be :state, :obs, :pseudo, or :shock"))
    end

    out = Base.filter(key -> dict[key] == index, collect(keys(dict)))
    if length(out) == 0
        error("Key corresponding to index $index not found for class :$class")
    elseif length(out) > 1
        error("Multiple keys corresponding to index $index found for class :$class")
    else
        return out[1]
    end
end

# Parse population mnemonic into 2 Nullable{Symbol}s from one
function parse_population_mnemonic(m::AbstractModel)
    mnemonic = get_setting(m, :population_mnemonic)
    if isnull(mnemonic)
        return [Nullable{Symbol}(), Nullable{Symbol}()]
    else
        return map(s -> Nullable(Symbol(s)), split(string(get(mnemonic)), DSGE_DATASERIES_DELIM))
    end
end

# From an augmented state space with anticipated policy shocks, get indices
# corresponding to pre-ZLB states, shocks, and observables
function inds_states_no_ant(m::AbstractModel)
    if n_anticipated_shocks(m) > 0
        ind_ant1 = m.endogenous_states[:rm_tl1]
        ind_antn = m.endogenous_states[Symbol("rm_tl$(n_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_states_augmented(m)]
    else
        return collect(1:n_states_augmented(m))
    end
end

function inds_shocks_no_ant(m::AbstractModel)
    if n_anticipated_shocks(m) > 0
        ind_ant1 = m.exogenous_shocks[:rm_shl1]
        ind_antn = m.exogenous_shocks[Symbol("rm_shl$(n_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_shocks_exogenous(m)]
    else
        return collect(1:n_shocks_exogenous(m))
    end
end

function inds_obs_no_ant(m::AbstractModel)
    if n_anticipated_shocks(m) > 0
        ind_ant1 = m.observables[:obs_nominalrate1]
        ind_antn = m.observables[Symbol("obs_nominalrate$(n_anticipated_shocks(m))")]
        return [1:(ind_ant1-1); (ind_antn+1):n_observables(m)]
    else
        return collect(1:n_observables(m))
    end
end

# Interface for I/O settings
spec(m::AbstractModel)         = m.spec
subspec(m::AbstractModel)      = m.subspec
saveroot(m::AbstractModel)     = get_setting(m, :saveroot)
dataroot(m::AbstractModel)     = get_setting(m, :dataroot)

# Interface for data
data_vintage(m::AbstractModel)    = get_setting(m, :data_vintage)
cond_vintage(m::AbstractModel)    = get_setting(m, :cond_vintage)
cond_id(m::AbstractModel)         = get_setting(m, :cond_id)
cond_full_names(m::AbstractModel) = get_setting(m, :cond_full_names)
cond_semi_names(m::AbstractModel) = get_setting(m, :cond_semi_names)
use_population_forecast(m::AbstractModel) = get_setting(m, :use_population_forecast)
hpfilter_population(m::AbstractModel)     = get_setting(m, :hpfilter_population)

# Interface for general computation settings
use_parallel_workers(m::AbstractModel)    = get_setting(m, :use_parallel_workers)

# Interface for estimation settings
reoptimize(m::AbstractModel)          = get_setting(m, :reoptimize)
calculate_hessian(m::AbstractModel) = get_setting(m, :calculate_hessian)
hessian_path(m::AbstractModel)      = get_setting(m, :hessian_path)
n_hessian_test_params(m::AbstractModel) = get_setting(m, :n_hessian_test_params)

# Interface for Metropolis-Hastings settings
n_mh_blocks(m::AbstractModel)      =  get_setting(m, :n_mh_blocks)
n_mh_simulations(m::AbstractModel) =  get_setting(m, :n_mh_simulations)
n_mh_burn(m::AbstractModel)        =  get_setting(m, :n_mh_burn)
mh_thin(m::AbstractModel)          =  get_setting(m, :mh_thin)

# Interface for forecast settings
date_forecast_start(m::AbstractModel)   = get_setting(m, :date_forecast_start)
forecast_block_size(m::AbstractModel)   = get_setting(m, :forecast_block_size)
forecast_start_block(m::AbstractModel)  = get_setting(m, :forecast_start_block)
forecast_input_file_overrides(m::AbstractModel) = get_setting(m, :forecast_input_file_overrides)
forecast_pseudoobservables(m::AbstractModel) = get_setting(m, :forecast_pseudoobservables)
forecast_uncertainty_override(m::AbstractModel) = get_setting(m, :forecast_uncertainty_override)
forecast_smoother(m::AbstractModel)     = get_setting(m, :forecast_smoother)
forecast_tdist_df_val(m::AbstractModel) = get_setting(m, :forecast_tdist_df_val)
forecast_tdist_shocks(m::AbstractModel) = get_setting(m, :forecast_tdist_shocks)
forecast_zlb_value(m::AbstractModel)    = get_setting(m, :forecast_zlb_value)
impulse_response_horizons(m::AbstractModel) = get_setting(m, :impulse_response_horizons)
n_shockdec_periods(m::AbstractModel)    = index_shockdec_end(m) - index_shockdec_start(m) + 1

function date_forecast_end(m::AbstractModel)
    date = date_forecast_start(m) + Dates.Month(3 * (forecast_horizons(m)-1))
    return Dates.lastdayofquarter(date)
end

function forecast_horizons(m::AbstractModel; cond_type::Symbol = :none)
    horizons = get_setting(m, :forecast_horizons)
    if cond_type == :none
        return horizons
    else
        return horizons - n_conditional_periods(m)
    end
end

function date_shockdec_start(m::AbstractModel)
    startdate = get_setting(m, :shockdec_startdate)
    if !isnull(startdate)
        return get(startdate)
    else
        return date_mainsample_start(m)
    end
end

function date_shockdec_end(m::AbstractModel)
    enddate =  get_setting(m, :shockdec_enddate)
    if !isnull(enddate)
        return get(enddate)
    else
        return date_forecast_end(m)
    end
end

"""
```
load_parameters_from_file(m::AbstractModel,path::String)
```
Returns a vector of parameters, read from a file, suitable for updating `m`.
"""
function load_parameters_from_file(m::AbstractModel, path::String)

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

    @assert length(x) == length(m.parameters)
    @assert eltype(x) == typeof(m.parameters[1].value)
    return x
end

"""
```
specify_mode!(m::AbstractModel, mode_file::String=""; verbose=:low)
```

Updates the values of `m.parameters` with the values from
`mode_file`. Sets `reoptimize` setting to `false`.

Usage: should be run before calling `estimate(m)`, e.g.:

    m = Model990()
    specify_mode!(m, modefile)
    estimate(m)
"""
function specify_mode!(m::AbstractModel, mode_file::String = ""; verbose=:low)

    m <= Setting(:reoptimize, false)

    if mode_file == ""
        mode_file = inpath(m, "user", "paramsmode.h5")
    end
    mode_file = normpath(mode_file)

    update!(m,load_parameters_from_file(m,mode_file))

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Loaded previous mode from $mode_file.")
    end

end

"""
```
specify_hessian(m::AbstractModel, path::String=""; verbose=:low)
```

Specify a Hessian matrix calculated at the posterior mode to use in the model estimation. If
no path is provided, will attempt to detect location.
"""
function specify_hessian(m::AbstractModel, path::String=""; verbose=:low)
    if isempty(path)
        path = inpath(m, "user", "hessian.h5")
    end

    if isfile(path) && splitext(path)[2] == ".h5"
        m <= Setting(:hessian_path, normpath(abspath(path)))
    else
        error("Invalid input Hessian file: $path",)
    end

    m <= Setting(:calculate_hessian, false)

    if VERBOSITY[verbose] >= VERBOSITY[:low]
        println("Specified hessian from $path.")
    end
end


#=
Build paths to where input/output/results data are stored.

Description:
Creates the proper directory structure for input and output files, treating the DSGE/save
    directory as the root of a savepath directory subtree. Specifically, the following
    structure is implemented:

    dataroot/

    savepathroot/
                 output_data/<spec>/<subspec>/log/
                 output_data/<spec>/<subspec>/<out_type>/raw/
                 output_data/<spec>/<subspec>/<out_type>/work/
                 output_data/<spec>/<subspec>/<out_type>/tables/
                 output_data/<spec>/<subspec>/<out_type>/figures/

Note: we refer to the savepathroot/output_data/<spec>/<subspec>/ directory as saveroot.
=#
# """
# ```
# logpath(model)
# ```
# Returns path to log file. Path built as
# ```
# <output root>/output_data/<spec>/<subspec>/log/log_<filestring>.log
# ```
# """
# function logpath(m::AbstractModel)
#     return savepath(m, "log", "log.log")
# end

strs = [:work, :raw, :tables, :figures, :log]
fns = [Symbol(x, "path") for x in strs]
for (str, fn) in zip(strs, fns)
    @eval begin
        # First eval function
        function $fn(m::AbstractModel,
                     out_type::String,
                     file_name::String = "",
                     filestring_addl::Vector{String}=Vector{String}())
            return savepath(m, out_type, $(string(str)), file_name, filestring_addl)
        end

        # Then, add docstring to it
        @doc $(
        """
        ```
        $fn(m::AbstractModel, out_type::String, file_name::String="", filestring_addl::Vector{String}=Vector{String}())
        ```

        Returns path to specific $str output file, creating containing directory as needed. If
        `file_name` not specified, creates and returns path to containing directory only. Path built
        as
        ```
        <output root>/output_data/<spec>/<subspec>/<out_type>/$str/<file_name>_<filestring>.<ext>
        ```
        """
        ) $fn
    end
end

# Not exposed to user. Actually create path and insert model string to file name.
function savepath(m::AbstractModel,
                  out_type::String,
                  sub_type::String,
                  file_name::String = "",
                  filestring_addl::Vector{String} = Vector{String}())
    # Containing directory
    dir = String(joinpath(saveroot(m), "output_data", spec(m), subspec(m), out_type, sub_type))

    if !isempty(file_name)
        base = filestring_base(m)
        return savepath(dir, file_name, base, filestring_addl)
    else
        return dir
    end
end

function savepath(dir::String,
                  file_name::String = "",
                  filestring_base::Vector{String} = Vector{String}(),
                  filestring_addl::Vector{String} = Vector{String}())
    if !isdir(dir)
        mkpath(dir)
    end

    if !isempty(file_name)
        (base, ext) = splitext(file_name)
        myfilestring = filestring(filestring_base, filestring_addl)
        file_name_detail = base * myfilestring * ext

        return joinpath(dir, file_name_detail)
    else
        return dir
    end
end


# Input data handled slightly differently, because it is not model-specific.
"""
```
inpath{T<:String}(m::AbstractModel, in_type::T, file_name::T="")
```

Returns path to specific input data file, creating containing directory as needed. If
`file_name` not specified, creates and returns path to containing directory only. Valid
`in_type` includes:

* `\"data\"`: recorded data
* `\"cond\"`: conditional data - nowcasts for the current forecast quarter, or related
* `\"user\"`: user-supplied data for starting parameter vector, hessian, or related
* `\"scenarios\"`: alternative scenarios

Path built as
```
<data root>/<in_type>/<file_name>
```
"""
function inpath(m::AbstractModel, in_type::String, file_name::String="")
    path = dataroot(m)
    # Normal cases.
    if in_type in ["data", "cond", "scenarios"]
        path = joinpath(path, in_type)
    # User-provided inputs. May treat this differently in the future.
    elseif in_type == "user"
        path = joinpath(path, "user")
    else
        error("Invalid in_type: ", in_type)
    end

    # Containing dir
    if !isdir(path)
        mkpath(path)
    end

    # If file_name provided, return full path
    if !isempty(file_name)
        path = joinpath(path, file_name)
    end

    return path
end

function filestring_base(m::AbstractModel)
    if !m.testing
        base = Vector{String}()
        for (skey, sval) in m.settings
            if sval.print
                push!(base, to_filestring(sval))
            end
        end
        return base
    else
        return ["test"]
    end
end

filestring(m::AbstractModel) = filestring(m, Vector{String}())
filestring(m::AbstractModel, d::String) = filestring(m, [String(d)])
function filestring(m::AbstractModel, d::Vector{String})
    base = filestring_base(m)
    return filestring(base, d)
end

function filestring(base::Vector{String}, d::Vector{String})
    filestrings = vcat(base, d)
    sort!(filestrings)
    return "_" * join(filestrings, "_")
end

function filestring(d::Vector{String})
    sort!(d)
    return "_" * join(d, "_")
end

"""
```
transform_to_model_space!{T<:AbstractFloat}(m::AbstractModel, values::Vector{T})
```

Transforms `values` from the real line to the model space, and assigns `values[i]` to
`m.parameters[i].value` for non-steady-state parameters. Recomputes the steady-state
paramter values.

### Arguments
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.
"""
function transform_to_model_space!{T<:AbstractFloat}(m::AbstractModel, values::Vector{T})
    new_values = transform_to_model_space(m.parameters, values)
    update!(m, new_values)
    steadystate!(m)
end

"""
```
update!{T<:AbstractFloat}(m::AbstractModel, values::Vector{T})
```

Update `m.parameters` with `values`, recomputing the steady-state parameter values.

### Arguments:
- `m`: the model object
- `values`: the new values to assign to non-steady-state parameters.
"""
function update!{T<:AbstractFloat}(m::AbstractModel, values::Vector{T})
    update!(m.parameters, values)
    steadystate!(m)
end

"""
```
rand(d::DegenerateMvNormal, m::AbstractModel; cc::AbstractFloat = 1.0)
```

Generate a draw from `d` with variance optionally scaled by `cc^2`.
"""
function rand(d::DegenerateMvNormal, m::AbstractModel; cc::AbstractFloat = 1.0)
    return d.μ + cc*d.σ*randn(m.rng, length(d))
end

"""
`rand_prior(m::AbstractModel; ndraws::Int = 100_000)`

Draw a random sample from the model's prior distribution.
"""
function rand_prior(m::AbstractModel; ndraws::Int = 100_000)
    T = typeof(m.parameters[1].value)
    npara = length(m.parameters)
    priorsim = Array{T}(ndraws, npara)

    for i in 1:ndraws
        priodraw = Array{T}(npara)

        # Parameter draws per particle
        for j in 1:length(m.parameters)

            priodraw[j] = if !m.parameters[j].fixed
                prio = rand(m.parameters[j].prior.value)

                # Resample until all prior draws are within the value bounds
                while !(m.parameters[j].valuebounds[1] < prio < m.parameters[j].valuebounds[2])
                    prio = rand(m.parameters[j].prior.value)
                end

                prio
            else
                m.parameters[j].value
            end
        end
        priorsim[i,:] = priodraw'
    end

    priorsim
end
