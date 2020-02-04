"""
```
DSGEVAR{T} <: AbstractVARModel{T}
```

"""
mutable struct DSGEVAR{T} <: AbstractVARModel{T}
    dsge::AbstractDSGEModel{T}
    observables::OrderedDict{Symbol,Int}
    shocks::OrderedDict{Symbol,Int}
    lags::Int
    spec::String
    subspec::String
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
    m = DSGEVAR{T}(dsge, OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), 0, spec, subspec)

    # Set λ for DSGEVAR estimation, by default, to 0.5
    m <= Setting(:λ, 0.5, true, "lambda", "weight on DSGE prior")

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
get_dsge(m::DSGEVAR)        = m.dsge
n_observables(m::DSGEVAR)   = length(m.observables)
n_shocks(m::DSGEVAR)        = length(m.shocks)

# Helper to check for valid VAR systems
function check_valid_dsgevar(m::DSGEVAR;
                             observables::Vector{Symbol} = collect(keys(m.observables)),
                             shocks::Vector{Symbol} = collect(keys(m.shocks)),
                             lags::Int = n_lags(m))
    # Check lags
    @assert lags >= 0 "Number of lags must be non-negative."

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
function update!(m::DSGEVAR; observables::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 shocks::Vector{Symbol} = Vector{Symbol}(undef, 0),
                 lags::Int = 0, check_valid::Bool = true)

    if check_valid
        check_valid_dsgevar(m; observables = observables, shocks = shocks, lags = lags)
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
                 lags::Int = 0, check_valid::Bool = true)

    m.dsge = dsge
    update!(m; observables = observables, shocks = shocks, lags = lags,
            check_valid = check_valid)

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
