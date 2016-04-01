abstract AbstractModel{T}

function Base.show{T<:AbstractModel}(io::IO, m::T)
    @printf io "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "%s\n" T
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

@inline function Base.setindex!{T<:Number}(m::AbstractModel, value::T, i::Integer)
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
setindex!{T<:AbstractParameter}(m::AbstractModel, param::T, i::Integer)
```

If `i`<length(m.parameters), overwrites m.parameters[i] with
param. Otherwise, overwrites m.steady_state[i-length(m.parameters).
"""
@inline function Base.setindex!{T<:AbstractParameter}(m::AbstractModel, param::T, i::Integer)
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

# Padding for number of anticipated policy shocks
n_anticipated_shocks_padding(m::AbstractModel) = get_setting(m, :n_anticipated_shocks_padding)

# First period we should incorporate zero bound expectations
# ZLB expectations should begin in 2008 Q4
function zlb_start_index!(m::AbstractModel, data::DataFrame)
    date = get_setting(m, :zlb_start_date)
    zlb_start = find(x->x==date, data[:date])[1]
    m <= Setting(:zlb_start_index, zlb_start, "Index of first period to incorporate zero bound expectation")

    zlb_start
end

zlb_start_index(m::AbstractModel) = get_setting(m, :zlb_start_index)

# Number of presample periods
n_presample_periods(m::AbstractModel) = get_setting(m, :n_presample_periods)

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

# Interface for I/O settings
spec(m::AbstractModel)         = m.spec
subspec(m::AbstractModel)      = m.subspec
saveroot(m::AbstractModel)     = get_setting(m, :saveroot)
dataroot(m::AbstractModel)     = get_setting(m, :dataroot)
data_vintage(m::AbstractModel) = get_setting(m, :data_vintage)

# Interface for general computation settings
use_parallel_workers(m::AbstractModel)    = get_setting(m, :use_parallel_workers)

# Interface for data step
use_population_forecast(m::AbstractModel) = get_setting(m, :use_population_forecast)

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


"""
```
load_parameters_from_file(m::AbstractModel,path::AbstractString)
```
Returns a vector of parameters, read from a file, suitable for updating `m`.
"""
function load_parameters_from_file(m::AbstractModel, path::AbstractString)

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
specify_mode!(m::AbstractModel, mode_file::AbstractString=""; verbose=:low)
```

Updates the values of `m.parameters` with the values from
`mode_file`. Sets `reoptimize` setting to `false`.

Usage: should be run before calling `estimate(m)`, e.g.:

    m = Model990()
    specify_mode!(m, modefile)
    estimate(m)
"""
function specify_mode!(m::AbstractModel, mode_file::AbstractString = ""; verbose=:low)

    m <= Setting(:reoptimize, false, "Optimize the posterior mode. If false, reads in mode from a file.")

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
specify_hessian(m::AbstractModel, path::AbstractString=""; verbose=:low)
```

Specify a Hessian matrix calculated at the posterior mode to use in the model estimation. If
no path is provided, will attempt to detect location.
"""
function specify_hessian(m::AbstractModel, path::AbstractString=""; verbose=:low)
    if isempty(path)
        path = inpath(m, "user", "hessian.h5")
    end

    if isfile(path) && splitext(path)[2] == ".h5"
        m <= Setting(:hessian_path, normpath(abspath(path)))
    else
        error("Invalid input Hessian file: $path",)
    end

    m <= Setting(:calculate_hessian, false, "Calculate the Hessian at the mode")

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
"""
```
logpath(model)
```
Returns path to log file. Path built as
```
<output root>/output_data/<spec>/<subspec>/log/log_<modelstring>.log
```
"""
function logpath(m::AbstractModel)
    return savepath(m, "log", "log.log")
end

strs = [:work, :raw, :tables, :figures]
fns = [symbol(x, "path") for x in strs]
for (str, fn) in zip(strs, fns)
    @eval begin
        # First eval function
        function $fn{T<:AbstractString}(m::AbstractModel,
                                             out_type::T,
                                             file_name::T="",
                                             model_string_addl::Vector{T}=Vector{T}())
            return savepath(m, out_type, $(string(str)), file_name, model_string_addl)
        end

        # Then, add docstring to it
        @doc $(
        """
        ```
        $fn{T<:AbstractString}(m::AbstractModel, out_type::T, file_name::T="")
        ```

        Returns path to specific $str output file, creating containing directory as needed. If
        `file_name` not specified, creates and returns path to containing directory only. Path built
        as
        ```
        <output root>/output_data/<spec>/<subspec>/<out_type>/$str/<file_name>_<modelstring>.<ext>
        ```
        """
        ) $fn
    end
end

# Not exposed to user. Actually create path and insert model string to file name.
function savepath{T<:AbstractString}(m::AbstractModel,
                                     out_type::T,
                                     sub_type::T,
                                     file_name::T="",
                                     model_string_addl::Vector{T}=Vector{T}())
    # Containing dir
    path = joinpath(saveroot(m), "output_data", spec(m), subspec(m), out_type, sub_type)
    if !isdir(path)
        mkpath(path)
    end

    # File with model string inserted
    if !isempty(file_name)
        if isempty(model_string_addl)
            model_string = modelstring(m)
        else
            model_string = modelstring(m, model_string_addl)
        end
        (base, ext) = splitext(file_name)
        file_name_detail = base * model_string * ext
        path = joinpath(path, file_name_detail)
    end

    return path
end

# Input data handled slightly differently, because it is not model-specific.
"""
```
inpath{T<:AbstractString}(m::AbstractModel, in_type::T, file_name::T="")
```

Returns path to specific input data file, creating containing directory as needed. If
`file_name` not specified, creates and returns path to containing directory only. Valid
`in_type` includes:

* `"data"`: recorded data
* `"cond"`: conditional data - nowcasts for the current forecast quarter, or related
* `"user"`: user-supplied data for starting parameter vector, hessian, or related

Path built as
```
<data root>/<in_type>/<file_name>
```
"""
function inpath{T<:AbstractString}(m::AbstractModel, in_type::T, file_name::T="")
    path = dataroot(m)
    # Normal cases.
    if in_type == "data" || in_type == "cond"
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

modelstring(m::AbstractModel) = modelstring(m, Vector{AbstractString}())
modelstring(m::AbstractModel, d::AbstractString) = modelstring(m, [d])
function modelstring{T<:AbstractString}(m::AbstractModel,
                                        d::Vector{T})
    if !m.testing
        filestrings = Vector{T}()
        for (skey, sval) in m.settings
            if sval.print
                push!(filestrings, to_filestring(sval))
            end
        end
        append!(filestrings, d)
        sort!(filestrings)
        return "_"*join(filestrings, "_")
    else
        return "_test"
    end
end

"""
```
transform_to_model_space!{T<:AbstractFloat}(m::AbstractModel, values::Vector{T})
```

Transforms `values` from the real line to the model space, and assigns `values[i]` to
`m.parameters[i].value` for non-steady-state parameters. Recomputes the steady-state
paramter values.

### Arguments
-`m`: the model object
-`values`: the new values to assign to non-steady-state parameters.
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
-`m`: the model object
-`values`: the new values to assign to non-steady-state parameters.
"""
function update!{T<:AbstractFloat}(m::AbstractModel, values::Vector{T})
    update!(m.parameters, values)
    steadystate!(m)
end

"""
```
rand{T<:AbstractFloat, U<:AbstractModel}(d::DegenerateMvNormal, m::U; cc::T = 1.0)
```

Generate a draw from d with variance optionally scaled by cc^2.
"""
function rand{T<:AbstractFloat, U<:AbstractModel}(d::DegenerateMvNormal, m::U; cc::T = 1.0)
    return d.μ + cc*d.σ*randn(m.rng, length(d))
end

