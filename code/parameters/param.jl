# We will have a prior file that looks like
#   α = Param(0.2,false,[-1,1],Beta(1,0.5))
#   β = Param(27, true, [-Inf,Inf], Gamma(1,0.5))
# We define Param to be a subtype of Number so that we can use numerical operation methods in
#   https://github.com/JuliaLang/julia/blob/master/base/promotion.jl
type Param <: Number
    value::Float64
    fixed::Bool
    bounds::(Float64, Float64)
    transf::Function
    prior::Distribution
    toreal::Function
    tomodel::Function

    function Param(value::Float64, fixed::Bool, bounds::(Float64, Float64), prior::Distribution, trtype::Int64, trbounds::(Float64, Float64); transf::Function = identity)
        if fixed
            prior = PointMass(value)
            trtype = 0
            trbounds = (value, value)
        end
        (a, b) = trbounds
        return new(value, fixed, bounds, transf, prior, toreal(trtype, a, b, 1.), tomodel(trtype, a, b, 1.))
    end
end

# Constructor for values given in getpara00_990.m as del = 0.025, for example
function Param(value::Float64)
    return Param(value, true, (value, value), PointMass(value), 0, (value, value))
end

# Methods so that arithmetic with parameters can be done tersely, like "θ.α + θ.β"
# Some cases will still require accessing the value field, e.g. "log(Θ.α.value)"
Base.convert(::Type{Float64}, α::Param) = α.value
Base.promote_rule(::Type{Float64}, ::Type{Param}) = Float64
Base.promote_rule(::Type{Param}, ::Type{Float64}) = Float64
Base.promote_rule(::Type{Int64}, ::Type{Param}) = Float64
Base.promote_rule(::Type{Param}, ::Type{Int64}) = Float64
+(α::Param, β::Param) = α.value + β.value
-(α::Param, β::Param) = α.value - β.value
*(α::Param, β::Param) = α.value * β.value
/(α::Param, β::Param) = α.value / β.value
^(α::Param, β::Param) = α.value ^ β.value
-(α::Param)           = -α.value



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
