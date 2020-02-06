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
