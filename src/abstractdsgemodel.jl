abstract AbstractDSGEModel{T<:AbstractFloat}

function Base.show{T<:AbstractDSGEModel}(io::IO, m::T)
    @printf io "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "%s\n" T
    @printf io "no. states:             %i\n" num_states(m)
    @printf io "no. anticipated shocks: %i\n" num_anticipated_shocks(m)
    @printf io "no. anticipated lags:   %i\n" num_anticipated_lags(m)
    @printf io "data vintage:           %s\n" data_vintage(m)
    @printf io "description:\n %s\n"          description(m)
end

@inline function Base.getindex(m::AbstractDSGEModel, i::Integer)
    if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]
    end
end

# need to define like this so we can disable bounds checking
@inline function Base.getindex(m::AbstractDSGEModel, k::Symbol)
    i = m.keys[k]
    @inbounds if i <= (j = length(m.parameters))
        return m.parameters[i]
    else
        return m.steady_state[i-j]
    end
end

@inline function Base.setindex!{T<:Number}(m::AbstractDSGEModel, value::T, i::Integer)
    if i <= (j = length(m.parameters))
        param = m.parameters[i]
        param.value = value
        if isa(param, ScaledParameter)
            param.scaledvalue = param.scaling(value)
        end
        return param
    else
        ssparam = m.steady_state[i-j]
        ssparam.value = value
        return ssparam
    end
end

"""
`setindex!{T<:AbstractParameter}(m::AbstractDSGEModel, param::T, i::Integer)`

If `i`<length(m.parameters), overwrites m.parameters[i] with
param. Otherwise, overwrites m.steady_state[i-length(m.parameters).
"""
@inline function Base.setindex!{T<:AbstractParameter}(m::AbstractDSGEModel, param::T, i::Integer)
    if i <= (j = length(m.parameters))
        m.parameters[i] = param
    else
        m.steady_state[i-j] = param
    end
    return param
end

Base.setindex!(m::AbstractDSGEModel, value, k::Symbol) = Base.setindex!(m, value, m.keys[k])


"""
(<=){T}(m::AbstractDSGEModel{T}, p::AbstractParameter{T})

Syntax for adding a parameter to a model: m <= parameter.
NOTE: If `p` is added to `m` and length(m.steady_state) > 0, `keys(m)` will not generate the index of `p` in `m.parameters`.
"""
function (<=){T}(m::AbstractDSGEModel{T}, p::AbstractParameter{T})

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

#=
"""
(<=){T}(m::AbstractDSGEModel{T}, ssp::SteadyStateParameter)

Add a new steady-state value to the model by appending `ssp` to the `m.steady_state` and adding `ssp.key` to `m.keys`.
"""
=#
function (<=){T}(m::AbstractDSGEModel{T}, ssp::SteadyStateParameter)

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

Distributions.logpdf(m::AbstractDSGEModel) = logpdf(m.parameters)
Distributions.pdf(m::AbstractDSGEModel) = exp(logpdf(m))

# Number of anticipated policy shocks
num_anticipated_shocks(m::AbstractDSGEModel) = get_setting(m, :num_anticipated_shocks)

# Padding for number of anticipated policy shocks
num_anticipated_shocks_padding(m::AbstractDSGEModel) = get_setting(m, :num_anticipated_shocks_padding)

# Number of periods back we should start incorporating zero bound expectations
# ZLB expectations should begin in 2008 Q4
num_anticipated_lags(m::AbstractDSGEModel) = get_setting(m, :num_anticipated_lags)

# TODO: This should be set when the data are read in
# Number of presample periods
num_presample_periods(m::AbstractDSGEModel) = get_setting(m, :num_presample_periods)

# Number of a few things that are useful 
num_states(m::AbstractDSGEModel)                 = length(m.endogenous_states)
num_states_augmented(m::AbstractDSGEModel)       = num_states(m) + length(m.endogenous_states_postgensys)
num_shocks_exogenous(m::AbstractDSGEModel)       = length(m.exogenous_shocks)
num_shocks_expectational(m::AbstractDSGEModel)   = length(m.expected_shocks)
num_equilibrium_conditions(m::AbstractDSGEModel) = length(m.equilibrium_conditions)
num_observables(m::AbstractDSGEModel)            = length(m.observables)
num_parameters(m::AbstractDSGEModel)             = length(m.parameters)
num_parameters_steady_state(m::AbstractDSGEModel)= length(m.steady_state)
num_parameters_free(m::AbstractDSGEModel)        = sum([!α.fixed for α in m.parameters])

# Interface for I/O settings
spec(m::AbstractDSGEModel)          = m.spec
subspec(m::AbstractDSGEModel)       = m.subspec 
modelpathroot(m::AbstractDSGEModel) = get_setting(m, :modelpathroot)
datapathroot(m::AbstractDSGEModel)  = get_setting(m, :datapathroot)
data_vintage(m::AbstractDSGEModel)  = get_setting(m, :data_vintage)

# Interface for general computation settings
use_parallel_workers(m::AbstractDSGEModel)    = get_setting(m, :use_parallel_workers)

# Interface for estimation settings
reoptimize(m::AbstractDSGEModel)          = get_setting(m, :reoptimize)
recalculate_hessian(m::AbstractDSGEModel) = get_setting(m, :recalculate_hessian)
max_hessian_free_params(m::AbstractDSGEModel) = get_setting(m, :max_hessian_free_params)

# Interface for Metropolis-Hastings settings
num_mh_blocks(m::AbstractDSGEModel)      =  get_setting(m, :num_mh_blocks)
num_mh_simulations(m::AbstractDSGEModel) =  get_setting(m, :num_mh_simulations) 
num_mh_burn(m::AbstractDSGEModel)        =  get_setting(m, :num_mh_burn)
mh_thinning_step(m::AbstractDSGEModel)   =  get_setting(m, :mh_thinning_step)

#=
"""
Build paths to where input/output/results data are stored.

Description:
Creates the proper directory structure for input and output files, treating the DSGE/save
    directory as the root of a savepath directory subtree. Specifically, the following
    structure is implemented:

    datapathroot/
                 
    savepathroot/
                 output_data/<spec>/<subspec>/log/
                 output_data/<spec>/<subspec>/<out_type>/raw/
                 output_data/<spec>/<subspec>/<out_type>/work/
                 output_data/<spec>/<subspec>/<out_type>/tables/
                 output_data/<spec>/<subspec>/<out_type>/figures/

Note: we refer to the savepathroot/output_data/<spec>/<subspec>/ directory as modelpathroot.
"""
=#
function logpath(m::AbstractDSGEModel)
    return modelpath(m, "log", "log.log")
end
function rawpath{T<:AbstractString}(m::AbstractDSGEModel, out_type::T, file_name::T)
    return modelpath(m, out_type, "raw", file_name)
end
function workpath{T<:AbstractString}(m::AbstractDSGEModel, out_type::T, file_name::T)
    return modelpath(m, out_type, "work", file_name)
end
function tablespath{T<:AbstractString}(m::AbstractDSGEModel, out_type::T, file_name::T)
    return modelpath(m, out_type, "tables", file_name)
end
function figurespath{T<:AbstractString}(m::AbstractDSGEModel, out_type::T, file_name::T)
    return modelpath(m, out_type, "figures", file_name)
end
    
function modelpath{T<:AbstractString}(m::AbstractDSGEModel, out_type::T, sub_type::T,
    file_name::T)

    # Containing dir
    path = joinpath(modelpathroot(m), "output_data", spec(m), subspec(m), out_type, sub_type)
    if !isdir(path) 
        mkpath(path) 
    end

    # File with model string inserted
    model_string = modelstring(m)
    (base, ext) = splitext(file_name)
    file_name_detail = base * "_" * model_string * ext
    path = joinpath(path, file_name_detail)

    return path
end

# Input data handled slightly differently, because it is not model-specific.
function inpath{T<:AbstractString}(m::AbstractDSGEModel, in_type::T, file_name::T)
    path = datapathroot(m)
    # Normal cases.
    if in_type == "data" || in_type == "cond"
        path = joinpath(path, in_type)
    end

    # User-provided inputs. May treat this differently in the future.
    if in_type == "user"
        path = joinpath(path, "user")
    end

    # Containing dir
    if !isdir(path)
        mkpath(path)
    end

    return joinpath(path, file_name)
end

function modelstring(m::AbstractDSGEModel)
    m.testing ? "_test" : join(values(m._filestrings),"_")
end

#=
doc"""
tomodel!{T<:AbstractFloat}(m::AbstractDSGEModel, values::Vector{T})

### Parameters:
-`m`: the model object
-`values`: the new values to assign to non-steady-state parameters.

### Description:
Transforms `values` from the real line to the model space, and assigns `values[i]` to `m.parameters[i].value` for non-steady-state parameters. Recomputes the steady-state paramter values.
"""
=#
function tomodel!{T<:AbstractFloat}(m::AbstractDSGEModel, values::Vector{T})
    new_values = tomodel(m.parameters, values)
    update!(m, new_values)
    return steadystate!(m)
end

#=
doc"""
update!{T<:AbstractFloat}(m::AbstractDSGEModel, values::Vector{T})

### Parameters:
-`m`: the model object
-`values`: the new values to assign to non-steady-state parameters.

### Description:
Update `m.parameters` with `values`, recomputing the steady-state parameter values.
"""
=#
function update!{T<:AbstractFloat}(m::AbstractDSGEModel, values::Vector{T})
    update!(m.parameters, values)
    return steadystate!(m) 
end

"""
rand{T<:AbstractFloat, U<:AbstractDSGEModel}(d::DegenerateMvNormal, m::U; cc::T = 1.0)

Generate a draw from d with variance optionally scaled by cc^2.
"""
function rand{T<:AbstractFloat, U<:AbstractDSGEModel}(d::DegenerateMvNormal, m::U; cc::T = 1.0)
    return d.μ + cc*d.σ*randn(m.rng, length(d))
end
    
