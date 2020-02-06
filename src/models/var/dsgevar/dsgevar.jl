"""
```
DSGEVAR{T} <: AbstractDSGEVARModel{T}
```
"""
mutable struct DSGEVAR{T} <: AbstractDSGEVARModel{T}
    dsge::AbstractDSGEModel{T}
    observables::OrderedDict{Symbol,Int}
    shocks::OrderedDict{Symbol,Int}
    lags::Int
    λ::T
    spec::String
    subspec::String
    testing::Bool
end

function Base.show(io::IO, m::DSGEVAR)
    @printf io "DSGE-VAR Model\n"
    @printf io "observables:      %i\n" n_observables(m)
    @printf io "data_vintage:     %s\n" data_vintage(m.dsge)
    @printf io "DSGE model:       %s\n" spec(get_dsge(m))
    @printf io "DSGE description: %s\n" description(get_dsge(m))
end

function DSGEVAR(dsge::AbstractDSGEModel{T}, shocks::Vector{Symbol}, subspec::String = "ss0";
                 custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                 copy_dsge::Bool = false, testing = false) where {T<:Real}

    # Initialize specs
    spec     = "dsgevar_" * ModelConstructors.spec(dsge)
    subspec  = subspec

    # Initialize empty DSGEVAR
    if copy_dsge
        dsge = deepcopy(dsge)
    end
    m = DSGEVAR{T}(dsge, OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), 0,
                   0., spec, subspec, testing)

    # Initialize shocks from DSGE
    update!(m; shocks = shocks, check_valid = false)

    # Initialize subspec
    init_subspec!(m)

    # Do checks of observables, shocks, and lags
    check_valid_dsgevar(m)

    return m
end

# Access and size functions
n_lags(m::DSGEVAR)          = m.lags
get_observables(m::DSGEVAR) = collect(keys(m.observables))
get_shocks(m::DSGEVAR)      = collect(keys(m.shocks))
get_λ(m::DSGEVAR)           = m.λ
get_dsge(m::DSGEVAR)        = m.dsge
n_observables(m::DSGEVAR)   = length(m.observables)
n_shocks(m::DSGEVAR)        = length(m.shocks)

# Helper to check for valid VAR systems
function check_valid_dsgevar(m::DSGEVAR{T};
                             observables::Vector{Symbol} = collect(keys(m.observables)),
                             shocks::Vector{Symbol} = collect(keys(m.shocks)),
                             lags::Int = n_lags(m), λ::T = m.λ) where {T<:Real}
    # Check lags
    @assert lags >= 0 "Number of lags must be non-negative."

    # Check λ is non-negative
    @assert λ >= 0 "The weight on the DSGE prior λ must be non-negative."

    # Check observables
    missing_obs = Vector{Symbol}(undef, 0)
    for k in observables
        if !(k in keys(m.dsge.observables)) && !(k in keys(m.dsge.pseudo_observables))
            push!(missing_obs, k)
        end
    end

    if !isempty(missing_obs)
        error("The following observables are not found in the underlying DSGE: \n" *
              join(string.(missing_obs), ", "))
    end

    # Check shocks
    missing_shocks = Vector{Symbol}(undef, 0)
    for k in shocks
        if !(k in keys(m.dsge.exogenous_shocks))
            push!(missing_shocks, k)
        end
    end

    if !isempty(missing_shocks)
        error("The following shocks are not found in the underlying DSGE: \n" *
              join(string.(missing_shocks), ", "))
    end
end

# Update VAR system, e.g. observables
function update!(m::DSGEVAR{T}; observables::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 shocks::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 lags::Int = 0, λ::T = m.λ, check_valid::Bool = true) where {T<:Real}

    if check_valid
        check_valid_dsgevar(m; observables = observables, shocks = shocks,
                            lags = lags, λ = λ)
    end

    if m.λ != λ
        m.λ = λ
    end
    if !isempty(observables)
        for k in keys(m.observables)
            delete!(m.observables, k)
        end
        for (i,k) in enumerate(observables)
            m.observables[k] = i
        end
    end
    if !isempty(shocks)
        for k in keys(m.shocks)
            delete!(m.shocks, k)
        end
        for (i,k) in enumerate(shocks)
            m.shocks[k] = i
        end
    end
    if lags > 0
        m.lags = lags
    end

    return m
end

function update!(m::DSGEVAR, dsge::AbstractDSGEModel;
                 observables::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 shocks::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 lags::Int = 0, λ::Float64 = m.λ, check_valid::Bool = true)

    m.dsge = dsge
    update!(m; observables = observables, shocks = shocks, lags = lags,
            λ = λ, check_valid = check_valid)

    return m
end

# Updating parameters
function update!(m::DSGEVAR{T}, values::Vector{T}) where {T<:Real}
    DSGE.update!(m.dsge, values)
    steadystate!(m.dsge)
end

function update!(m::DSGEVAR{T}, values::ParameterVector{T}) where {T<:Real}
    DSGE.update!(m.dsge, values)
    steadystate!(m.dsge)
end

# Overload <= so that it applies to the underlying DSGE object
function (<=)(m::DSGEVAR{T}, p::ModelConstructors.AbstractParameter{T}) where {T}
    m.dsge <= p
end

function (<=)(m::DSGEVAR{T}, p::Union{ModelConstructors.SteadyStateParameter, ModelConstructors.SteadyStateParameterArray}) where {T}
    m.dsge <= p
end

function (<=)(m::DSGEVAR{T}, p::ModelConstructors.SteadyStateParameterGrid) where {T}
    m.dsge <= p
end

function (<=)(m::DSGEVAR, s::Setting)
    m.dsge <= s
end

# Setting access
function get_setting(m::DSGEVAR, k::Symbol)
    return ModelConstructors.get_setting(m.dsge, k)
end

# Prior
function prior(m::DSGEVAR)
    return ModelConstructors.prior(m.dsge)
end

# Not exposed to user. Actually create path and insert model string to file name.
function savepath(m::DSGEVAR,
                  out_type::String,
                  sub_type::String,
                  file_name::String = "",
                  filestring_addl::Vector{String} = Vector{String}())
    # Containing directory
    dir = String(joinpath(saveroot(m), "output_data", spec(m.dsge), subspec(m.dsge),
                          spec(m), subspec(m), out_type, sub_type))

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

function filestring_base(m::DSGEVAR)
    if !m.testing
        base = Vector{String}()
        for (skey, sval) in m.dsge.settings
            if sval.print
                push!(base, ModelConstructors.to_filestring(sval))
            end
        end
        return base
    else
        return ["test"]
    end
end
