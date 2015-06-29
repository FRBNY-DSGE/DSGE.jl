using Distributions: Distribution
import Base: convert, promote_rule, log, exp

using DSGE.DistributionsExt: PointMass

# We define Param to be a subtype of Number so that we can use numerical operation methods in
#   https://github.com/JuliaLang/julia/blob/master/base/promotion.jl

# Some parameter values, especially rates, are given in one form but in practice we would like to
# use them in another form. For example, pistar = 0.5000 refers to an inflation rate of 0.5%; its
# transf function x -> 1 + x/100 specifies the form in which pistar is used in calculations. We
# should ALWAYS use the transformed parameter value (tval) in calculations.

type Param <: Number
    value::Float64
    transf::Function
    tval::Float64
    fixed::Bool
    bounds::(Float64, Float64)
    prior::Distribution
    trtype::Int64
    trbounds::(Float64, Float64)
    description::String

    function Param(value::Float64, fixed::Bool, bounds::(Float64, Float64), prior::Distribution,
                   trtype::Int64, trbounds::(Float64, Float64); transf::Function = identity,
                   description::String = "")
        if fixed
            prior = PointMass(value)
            trtype = 0
            trbounds = (value, value)
        end
        if trtype != 0 && trtype != 1 && trtype != 2
            error("trtype must be 0, 1, or 2")
        end
        (a, b) = trbounds
        return new(value, transf, transf(value), fixed, bounds, prior, trtype, trbounds, description)
    end
end

# Constructor for values given in getpara00_990.m as del = 0.025, for example
function Param(value::Float64)
    return Param(value, true, (value, value), PointMass(value), 0, (value, value))
end



# Methods so that arithmetic with parameters can be done tersely, like "θ.α + θ.β"
# Note there are still instances where we must refer to α.tval, e.g. pdf(α.prior, α.val)
Base.convert(::Type{Float64}, α::Param) = α.tval
Base.promote_rule{T<:FloatingPoint}(::Type{Param}, ::Type{T}) = Float64
Base.promote_rule{T<:Integer}(::Type{Param}, ::Type{T}) = Float64

for op in [:+, :-, :*, :/, :^]
    @eval ($op)(α::Param, β::Param) = ($op)(α.tval, β.tval)
end

for f in [:-, :log, :exp]
    @eval ($f)(α::Param) = $(f)(α.tval)
end

# Transforms variables from model to max (invtrans.m)
function toreal(α::Param)
    (a, b) = α.trbounds
    c = 1.

    if α.trtype == 0
       return α.tval
    elseif α.trtype == 1
        cx = 2 * (α.tval - (a+b)/2) / (b-a);
        return (1/c) * cx / sqrt(1 - cx^2);
    elseif α.trtype == 2
        return b + (1/c) * log(α.tval-a)
    end
end

# Transforms variables from max to model (trans.m)
function tomodel(α::Param)
    (a, b) = α.trbounds
    c = 1.

    if α.trtype == 0
        return α.tval
    elseif α.trtype == 1
        return (a+b)/2 + (b-a)/2 * c * α.tval / sqrt(1 + c^2 * α.tval^2)
    elseif α.trtype == 2
        return a + exp(c * (α.tval-b))
    end
end
