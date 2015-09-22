import Base: convert, promote_rule, +, -, *, /, ^, log, exp
using Compat

"""
## Params

We define Param to be a subtype of Number so we can use numerical operation methods in
  https://github.com/JuliaLang/julia/blob/master/base/promotion.jl

Some parameter values, especially rates, are given in one form but in practice we would
like to use them in another form. For example, pistar = 0.5000 refers to an inflation rate
of 0.5%; its scalefunction x -> 1 + x/100 specifies the form in which pistar is used in
calculations. We use the scaled parameter value in calculations, but the original value in
estimation.
"""
:param_type

typealias Interval{T} @compat Tuple{T, T}

# Compat should take care of this
## if VERSION < v"0.4.0-rc+10"
##     typealias AbstractString String
## end

type Param <: Number
    value::Float64
    scalefunction::Function
    scaledvalue::Float64
    fixed::Bool
    bounds::Interval{Float64}
    priordist::Distribution
    transformtype::Int64
    transformbounds::Interval{Float64}
    description::AbstractString
    texLabel::AbstractString
    
    function Param{T<:AbstractString}(value::Float64, fixed::Bool, bounds::Interval{Float64},
                   priordist::Distribution, transformtype::Int64,
                   transformbounds::Interval{Float64}; scalefunction::Function = identity,
                   description::T = "",texLabel::T="")
        if fixed
            priordist = PointMass(value)
            transformtype = 0
        end
        if transformtype != 0 && transformtype != 1 && transformtype != 2
            error("transformtype must be 0, 1, or 2")
        end
        (a, b) = transformbounds
        return new(value, scalefunction, scalefunction(value), fixed, bounds, priordist,
                   transformtype, transformbounds, description,texLabel)
    end
end

# Constructor for values given in getpara00_990.m as del = 0.025, for example
function Param(value::Float64)
    return Param(value, true, (value, value), PointMass(value), 0, (value, value))
end

# Update a vector of parameters with a vector of new values.
# TODO the user should never be able to update the parameters without also calling
# steadystate!(model). See estimate:32,33. However, if we want to operate on vectors of
# parameters in the file parameters.jl, then we don't have access to the model object and
# the steadystate! function.
function update!{T<:AbstractFloat}(parameters::Vector{Param}, newvalues::Vector{T})
    @assert length(newvalues) == length(parameters)
    map(update!, parameters, newvalues)
    return parameters
end

# Update a Param's value and scaledvalue if it is not fixed
function update!(θ::Param, newvalue::Float64)
    if !θ.fixed
        θ.value = newvalue
        θ.scaledvalue = θ.scalefunction(newvalue)
    end
    return θ
end


# Methods so that arithmetic with parameters can be done tersely, like "θ.α + θ.β"
# Note there are still cases where we must refer to α.scaledvalue, e.g. pdf(θ.priordist, θ.val)
Base.convert{T<:AbstractFloat}(::Type{T}, α::Param) = α.scaledvalue
Base.promote_rule{T<:AbstractFloat}(::Type{Param}, ::Type{T}) = Float64
Base.promote_rule{T<:Integer}(::Type{Param}, ::Type{T}) = Float64

for op in [:+, :-, :*, :/, :^]
    @eval ($op)(α::Param, β::Param) = ($op)(α.scaledvalue, β.scaledvalue)
end

for f in [:-, :log, :exp]
    @eval ($f)(α::Param) = $(f)(α.scaledvalue)
end

# Returns a vector of parameter values transformed to lie on the real line.
function toreal(parameters::Vector{Param})
    return [toreal(θ) for θ in parameters]
end

# Transforms variables from model to max (invtrans.m)
function toreal(θ::Param)
    (a, b) = θ.transformbounds
    c = 1.0

    if θ.transformtype == 0
       return θ.value
    elseif θ.transformtype == 1
        cx = 2 * (θ.value - (a+b)/2) / (b-a)
        return (1/c) * cx / sqrt(1 - cx^2)
    elseif θ.transformtype == 2
        return b + (1/c) * log(θ.value-a)
    else
        error("Invalid transform type $θ.transformtype")
    end
end

function tomodel{T<:AbstractFloat}(values::Vector{T}, parameters::Vector{Param})
    return map(tomodel, values, parameters)
end

# Transforms variables from max to model (trans.m)
function tomodel{T<:AbstractFloat}(value::T, θ::Param)
    (a, b) = θ.transformbounds
    c = 1.0

    if θ.transformtype == 0
        return value
    elseif θ.transformtype == 1
        return (a+b)/2 + (b-a)/2*c*value/sqrt(1 + c^2 * value^2)
    elseif θ.transformtype == 2
        return a + exp(c * (value-b))
    else
        error("Invalid transform type $θ.transformtype")
    end
end

# Given a vector of parameter values on the real line, map them to the model space and
# update model.parameters field.
function tomodel!{T<:AbstractFloat}(values::Vector{T}, parameters::Vector{Param})
    newvalues = map(tomodel, values, parameters)
    return update!(parameters, newvalues)
end
