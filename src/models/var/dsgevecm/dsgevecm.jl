"""
```
DSGEVECM{T} <: AbstractDSGEVECMModel{T}
```
implements a simple interface for combining a given DSGE model
with a VECM (a VAR with cointegrating terms) to create a DSGE-VECM.
Confer with Del Negro and Schorfheide (2006)
and/or Del Negro, Schorfheide, Smets, and Wouters (2007) for details about DSGE-VECMs.
We recommend the first paper as an initial introduction to DSGE-VECMs.

The recommended constructor requires the user to provide
(1) an `AbstractDSGEModel` object, (2) which structural shocks from
the DSGE to use, and (3) the subspec (optional, defaults to "ss0").
If the subspec "ss0" is used, then the result is a `DSGEVECM`
whose VECM component is "empty" in that the
observables, cointegrating relationships, lags, and λ weight are not specified.
The reason why this constructor requires the user to specify which
structural shocks of DSGE to use is that this information is
DSGE-specific rather than information about the VECM.

However, we can also construct a `DSGEVECM` without having to
specify the structural shocks when calling the constructor,
although we still need to give an instance of an `AbstractDSGEModel`.

### Example
The code below instantiates an empty `DSGEVECM`
with `AnSchorfheide` as the underlying DSGE and then
calls `update!` on the empty DSGE-VECM
to add information about the desired DSGE-VECM spec.

```jldoctest; output = false
dsgevecm = DSGEVECM(AnSchorfheide())
DSGE.update!(dsgevecm, shocks = [:rm_sh, :z_sh, :g_sh],
    observables = [:obs_gdp, :obs_cpi, :obs_nominalrate],
    λ = 1., lags = 4)

# output

DSGE-VECM Model
observables:      3
data_vintage:     200310
DSGE model:       an_schorfheide
DSGE description: Julia implementation of model defined in 'Bayesian Estimation of DSGE Models' by Sungbae An and Frank Schorfheide: AnSchorfheide, ss0
```

### Fields

#### DSGE object
* `dsge::AbstractDSGEModel{T}`: underlying DSGE model object

#### DSGE-VECM Information
* `observables::OrderedDict{Symbol,Int}`: dictionary mapping observables
    of the VECM to their index in the matrices representing the DSGE-VECM.
* `cointegrating::OrderedDict{Symbol,Int}`: dictionary mapping cointegrating relationships
    of the VECM to their index in the matrices representing the DSGE-VECM.
    When creating the state space representation of a DSGE-VECM, the
    cointegrating relationships will come after the observables. Accordingly,
    the first `cointegrating` index will be after the last observable.
* `shocks::OrderedDict{Symbol,Int}`: dictionary mapping structural
    shocks in the DSGE to their index in the matrices representing the DSGE-VECM
* `lags::Int`: number of lags in the VECM
* `λ::T`: weight on the DSGE prior

#### Auxiliary Information
* `spec::String`: concatenates `dsgevar` with the spec of the DSGE, e.g.
    for `AnSchorfheide`, we have `dsgevar_an_schorfheide`.
* `subspec::String`: specifies the model subspecification.
    Cached here for filepath computation.
* `testing::Bool`: indicates whether the model is in testing mode.
    Currently, this setting has no uses in practice
"""
mutable struct DSGEVECM{T} <: AbstractDSGEVECMModel{T}
    dsge::AbstractDSGEModel{T}
    observables::OrderedDict{Symbol,Int}
    cointegrating::OrderedDict{Symbol,Int}
    shocks::OrderedDict{Symbol,Int}
    lags::Int
    λ::T
    spec::String
    subspec::String
    testing::Bool
end

function Base.show(io::IO, m::DSGEVECM)
    @printf io "DSGE-VECM Model\n"
    @printf io "observables:      %i\n" n_observables(m)
    @printf io "data_vintage:     %s\n" data_vintage(m.dsge)
    @printf io "DSGE model:       %s\n" spec(get_dsge(m))
    @printf io "DSGE description: %s\n" description(get_dsge(m))
end

function DSGEVECM(dsge::AbstractDSGEModel{T}, shocks::Vector{Symbol}, subspec::String = "ss0";
                 custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                 copy_dsge::Bool = false, testing = false) where {T<:Real}

    # Initialize specs
    spec     = "dsgevecm_" * ModelConstructors.spec(dsge)
    subspec  = subspec

    # Initialize empty DSGEVECM
    if copy_dsge
        dsge = deepcopy(dsge)
    end
    m = DSGEVECM{T}(dsge, OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),
                    OrderedDict{Symbol,Int}(), 0,
                    0., spec, subspec, testing)

    # Initialize shocks from DSGE
    update!(m; shocks = shocks, check_valid = false)

    # Initialize subspec
    init_subspec!(m)

    # Do checks of observables, shocks, and lags
    check_valid_dsgevecm(m)

    return m
end

# Empty constructor
function DSGEVECM(dsge::AbstractDSGEModel{T}, subspec::String = "ss0";
                 custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                 copy_dsge::Bool = false, testing = false) where {T<:Real}

    # Initialize specs
    spec     = "dsgevar_" * ModelConstructors.spec(dsge)
    subspec  = subspec

    # Initialize empty DSGEVECM
    if copy_dsge
        dsge = deepcopy(dsge)
    end
    m = DSGEVECM{T}(dsge, OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),
                    OrderedDict{Symbol,Int}(), 0,
                    0., spec, subspec, testing)

    # Initialize shocks from DSGE
    update!(m; check_valid = false)

    # Initialize subspec
    init_subspec!(m)

    # Do checks of observables, shocks, and lags
    check_valid_dsgevecm(m)

    return m
end
