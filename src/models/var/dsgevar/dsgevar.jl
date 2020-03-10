"""
```
DSGEVAR{T} <: AbstractDSGEVARModel{T}
```
implements a simple interface for combining a given DSGE model
with a VAR to create a DSGE-VAR. Confer with Del Negro and Schorfheide (2004),
Del Negro and Schorfheide (2006), and Del Negro and Schorfheide (2009) for
details about DSGE-VARs.

We can easily construct an "empty" DSGE-VAR with `AnSchorfheide` as
the underlying DSGE and then call `update!` on the empty DSGE-VAR
to add information about the desired DSGE-VAR spec.


```jldoctest; output = false
dsgevar = DSGEVAR(AnSchorfheide())
DSGE.update!(dsgevar, shocks = [:rm_sh, :z_sh, :g_sh],
    observables = [:obs_gdp, :obs_cpi, :obs_nominalrate],
    λ = 1., lags = 4)

# output

DSGE-VAR Model
observables:      3
data_vintage:     200310
DSGE model:       an_schorfheide
DSGE description: Julia implementation of model defined in 'Bayesian Estimation of DSGE Models' by Sungbae An and Frank Schorfheide: AnSchorfheide, ss0
```

When extensively using DSGE-VARs, we recommend defining your own subspecs in
`subspecs.jl` because it simplifies the process of saving, reading, and analyzing
output from estimating and calculating impulse responses for DSGE-VARs.

### Fields

#### DSGE object
* `dsge::AbstractDSGEModel{T}`: underlying DSGE model object

#### DSGE-VAR Information
* `observables::OrderedDict{Symbol,Int}`: dictionary mapping observables
    of the VAR to their index in the matrices representing the DSGE-VAR
* `shocks::OrderedDict{Symbol,Int}`: dictionary mapping structural
    shocks in the DSGE to their index in the matrices representing the DSGE-VAR
* `lags::Int`: number of lags in the VAR
* `λ::T`: weight on the DSGE prior

#### Auxiliary Information
* `spec::String`: concatenates `dsgevar` with the spec of the DSGE, e.g.
    for `AnSchorfheide`, we have `dsgevar_an_schorfheide`.
* `subspec::String`: specifies the model subspecification.
    Cached here for filepath computation.
* `testing::Bool`: indicates whether the model is in testing mode.
    Currently, this setting has no uses in practice
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

# Empty constructor
function DSGEVAR(dsge::AbstractDSGEModel{T}, subspec::String = "ss0";
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
    update!(m; check_valid = false)

    # Initialize subspec
    init_subspec!(m)

    # Do checks of observables, shocks, and lags
    check_valid_dsgevar(m)

    return m
end
