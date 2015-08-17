abstract AbstractModel

# We define Param to be a subtype of Number so we can use numerical operation methods in
#   https://github.com/JuliaLang/julia/blob/master/base/promotion.jl

# Some parameter values, especially rates, are given in one form but in practice we would
# like to use them in another form. For example, pistar = 0.5000 refers to an inflation rate
# of 0.5%; its scalefunction x -> 1 + x/100 specifies the form in which pistar is used in
# calculations. We use the scaled parameter value in calculations, but the original value in
# estimation.

typealias Interval{T} @compat Tuple{T, T}

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
    c = 1.0

    if α.transformtype == 0
       return α.value
    elseif α.transformtype == 1
        cx = 2 * (α.value - (a+b)/2) / (b-a)
        return (1/c) * cx / sqrt(1 - cx^2)
    elseif α.transformtype == 2
        return b + (1/c) * log(α.value-a)
    else
        error("Invalid transform type $α.transformtype")
    end
end

# Transforms variables from max to model (trans.m)
function tomodel{T<:FloatingPoint}(value::T, α::Param)
    (a, b) = α.transformbounds
    c = 1.0

    if α.transformtype == 0
        return value
    elseif α.transformtype == 1
        return (a+b)/2 + (b-a)/2*c*value/sqrt(1 + c^2 * value^2)
    elseif α.transformtype == 2
        return a + exp(c * (value-b))
    else
        error("Invalid transform type $α.transformtype")
    end
end

# The abstract Parameters type is the supertype of all model-specific ParametersXXX types.
# All concrete types have both Param (parameters) and Float64 (steady-state values) fields.
# See Parameters990 for an example.
abstract Parameters

# Implement the iterator protocol for the Parameters type
# This will iterate over all Param fields (not steady-state values)
Base.start(Θ::Parameters) = 1

function Base.next(Θ::Parameters, state::Int)
    α = getfield(Θ, state)
    state += 1
    while !done(Θ, state) && !isa(getfield(Θ, state), Param)
        state += 1
    end
    return α, state
end

Base.done(Θ::Parameters, state::Int) = state == length(names(Θ))+1

# Length of a Parameters object is the number of Param fields
Base.length(Θ::Parameters) = count(field -> isa(getfield(Θ, field), Param), names(Θ))

function update!{T<:FloatingPoint}(Θ::Parameters, newvalues::Vector{T})
    @assert length(newvalues) == length(Θ)
    for (α, newvalue) in zip(Θ, newvalues)
        update!(α, newvalue)
    end
    return steadystate!(Θ)
end

# Returns a vector of parameter values transformed to lie on the real line
function toreal(Θ::Parameters)
    return [toreal(α) for α in Θ]
end

# Given a vector of parameter values on the real line, maps them to the model space and updates Θ
function tomodel!{T<:FloatingPoint}(values::Vector{T}, Θ::Parameters)
    newvalues = [tomodel(value, α) for (value, α) in zip(values, Θ)]
    return update!(Θ, newvalues)
end

# Calculate (log of) joint density of Θ
function prior(Θ::Parameters)
    sum = 0.0
    for φ in Θ
        curr = logpdf(φ.priordist, φ.value)
        sum += curr
    end
    return sum
end

# A type that bundles together the model indices dictionaries from
#   models/m$(spec)/modelinds.jl
type ModelInds
    endostates::Dict{String, Int64}
    exoshocks::Dict{String, Int64}
    expshocks::Dict{String, Int64}
    eqconds::Dict{String, Int64}
    endostates_postgensys::Dict{String, Int64}
    observables::Dict{String, Int64}
end

# Given an array of names, return a dictionary mapping names to indices
# When optional field `start` provided, first index is start+1
function makedict{T<:String}(names::Vector{T}; start::Int = 0)
    return [names[i] => start+i for i = 1:length(names)]
end
