using DistributionsExt.PointMass

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
    toreal::Function
    tomodel::Function
    description::String

    function Param(value::Float64, fixed::Bool, bounds::(Float64, Float64), prior::Distribution,
                   trtype::Int64, trbounds::(Float64, Float64); transf::Function = identity,
                   description::String = "")
        if fixed
            prior = PointMass(value)
            trtype = 0
            trbounds = (value, value)
        end
        (a, b) = trbounds
        return new(value, transf, transf(value), fixed, bounds, prior, toreal(trtype, a, b, 1.),
                   tomodel(trtype, a, b, 1.), description)
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
+(α::Param, β::Param) = α.tval + β.tval
-(α::Param, β::Param) = α.tval - β.tval
*(α::Param, β::Param) = α.tval * β.tval
/(α::Param, β::Param) = α.tval / β.tval
^(α::Param, β::Param) = α.tval ^ β.tval
-(α::Param)           = -α.tval
Base.log(α::Param)    = log(α.tval)
Base.exp(α::Param)    = exp(α.tval)



# Transforms variables from model to max (invtrans.m)
# Returns an anonymous function (Float64 -> Float64)
function toreal(trtype::Int64, a::Float64, b::Float64, c::Float64)
    if trtype == 0
       return identity
    elseif trtype == 1
        return function (x)
            cx = 2 * (x - (a+b)/2) / (b-a);
            return (1/c) * cx / sqrt(1 - cx^2);
        end
    elseif trtype == 2
        return x -> b + (1/c) * log(x-a)
    else
        error("trtype must be 0, 1, or 2")
    end
end

# Transforms variables from max to model (trans.m)
function tomodel(trtype::Int64, a::Float64, b::Float64, c::Float64)
    if trtype == 0
        return identity
    elseif trtype == 1
        return x -> (a+b)/2 + (b-a)/2 * c * x / sqrt(1 + c^2 * x^2)
    elseif trtype == 2
        return x -> a + exp(c * (x-b))
    else
        error("trtype must be 0, 1, or 2")
    end
end
