abstract AbstractDSGEModel

Base.getindex(m::AbstractDSGEModel, i::Integer) = m.par[i]
Base.getindex(m::AbstractDSGEModel, keyword::Symbol) = m.par[m.parkeys[keyword]]

Base.setindex!(m::AbstractDSGEModel, value, i::Integer) = setindex!(m.par, value, i)
Base.setindex!(m::AbstractDSGEModel, value, keyword::Symbol) = setindex!(m.par, value, m.parkeys[keyword])

function Base.show{T<:AbstractDSGEModel}(io::IO, m::T)
    @printf io "%s" "Dynamic Stochastic General Equilibrium Model\n"
    @printf io "%s %s\n" "model: " T
    @printf io "%s             %i\n" "# states:" states(m)
    @printf io "%s %i\n" "# anticipated shocks:" anticipated_shocks(m)
    @printf io "%s   %i\n" "# anticipated lags:" anticipated_lags(m)
    @printf io "%s\n %s\n" "description:" description(m)
end

# Number of anticipated policy shocks
anticipated_shocks(m::AbstractDSGEModel) = m.anticipated_shocks

# Padding for nant
anticipated_shocks_padding(m::AbstractDSGEModel) = m.anticipated_shocks_padding

# Number of periods back we should start incorporating zero bound expectations
# ZLB expectations should begin in 2008 Q4
anticipated_lags(m::AbstractDSGEModel) = m.anticipated_lags

# TODO: This should be set when the data are read in
# Number of presample periods
presample_periods(m::AbstractDSGEModel) = m.presample_periods

# Number of a few things that are useful apparently
states(m::AbstractDSGEModel)      = length(m.endogenous_states)
augmented_states(m::AbstractDSGEModel)  = states(m) + length(m.endogenous_states_postgensys)
exogenous_shocks(m::AbstractDSGEModel)   = length(m.exogenous_shocks)
expected_shocks(m::AbstractDSGEModel)   = length(m.expected_shocks)
equilibrium_conditions(m::AbstractDSGEModel)     = length(m.equilibrium_conditions)
observables(m::AbstractDSGEModel) = length(m.observables)

# We define Param to be a subtype of Number so we can use numerical operation methods in
#   https://github.com/JuliaLang/julia/blob/master/base/promotion.jl

# Some parameter values, especially rates, are given in one form but in practice we would
# like to use them in another form. For example, pistar = 0.5000 refers to an inflation rate
# of 0.5%; its scalefunction x -> 1 + x/100 specifies the form in which pistar is used in
# calculations. We use the scaled parameter value in calculations, but the original value in
# estimation.

typealias Interval{T} @compat Tuple{T, T}

# The abstract Parameters type is the supertype of all model-specific ParametersXXX types.
# All concrete types have both Param (parameters) and Float64 (steady-state values) fields.
# See Parameters990 for an example.
abstract Parameters

# of states:             66
# of anticipated shocks: 6
# of anticipated lags:   24

type Param <: Number
    value::Float64
    scalefunction::Function
    scaledvalue::Float64
    fixed::Bool
    bounds::Interval{Float64}
    priordist::Distribution
    transformtype::Int64
    transformbounds::Interval{Float64}
    description::String

    function Param(value::Float64, fixed::Bool, bounds::Interval{Float64},
                   priordist::Distribution, transformtype::Int64,
                   transformbounds::Interval{Float64}; scalefunction::Function = identity,
                   description::String = "")
        if fixed
            priordist = PointMass(value)
            transformtype = 0
        end
        if transformtype != 0 && transformtype != 1 && transformtype != 2
            error("transformtype must be 0, 1, or 2")
        end
        (a, b) = transformbounds
        return new(value, scalefunction, scalefunction(value), fixed, bounds, priordist,
                   transformtype, transformbounds, description)
    end
end

# Constructor for values given in getpara00_990.m as del = 0.025, for example
function Param(value::Float64)
    return Param(value, true, (value, value), PointMass(value), 0, (value, value))
end

# Update a Param's value and scaledvalue if it is not fixed
function update!(α::Param, newvalue::Float64)
    if !α.fixed
        α.value = newvalue
        α.scaledvalue = α.scalefunction(newvalue)
    end
    return α
end

function update!{T<:FloatingPoint}(Θ::Parameters, newvalues::Vector{T})
    @assert length(newvalues) == length(Θ)
    for (α, newvalue) in zip(Θ, newvalues)
        update!(α, newvalue)
    end
    return steadystate!(Θ)
end

# Methods so that arithmetic with parameters can be done tersely, like "θ.α + θ.β"
# Note there are still cases where we must refer to α.scaledvalue, e.g. pdf(α.priordist, α.val)
Base.convert{T<:FloatingPoint}(::Type{T}, α::Param) = α.scaledvalue
Base.promote_rule{T<:FloatingPoint}(::Type{Param}, ::Type{T}) = Float64
Base.promote_rule{T<:Integer}(::Type{Param}, ::Type{T}) = Float64

for op in [:+, :-, :*, :/, :^]
    @eval ($op)(α::Param, β::Param) = ($op)(α.scaledvalue, β.scaledvalue)
end

for f in [:-, :log, :exp]
    @eval ($f)(α::Param) = $(f)(α.scaledvalue)
end

# Transforms variables from model to max (invtrans.m)
function toreal(α::Param)
    (a, b) = α.transformbounds
    c = 1.

    if α.transformtype == 0
       return α.scaledvalue
    elseif α.transformtype == 1
        cx = 2 * (α - (a+b)/2) / (b-a)
        return (1/c) * cx / sqrt(1 - cx^2)
    elseif α.transformtype == 2
        return b + (1/c) * log(α-a)
    end
end

# Transforms variables from max to model (trans.m)
function tomodel(α::Param)
    (a, b) = α.transformbounds
    c = 1.

    if α.transformtype == 0
        return α.scaledvalue
    elseif α.transformtype == 1
        return (a+b)/2 + (b-a)/2*c*α/sqrt(1 + c^2 * α.scaledvalue^2)
    elseif α.transformtype == 2
        return a + exp(c * (α-b))
    end
end
