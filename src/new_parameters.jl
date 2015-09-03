typealias Interval{T} @compat Tuple{T,T}

# define all the kinds of transformations we make
abstract Transform
immutable Untransformed <: Transform end
immutable SquareRoot    <: Transform end
immutable Exponential   <: Transform end

Base.show(io::IO, t::Untransformed) = @printf io "x -> x\n"
Base.show(io::IO, t::SquareRoot)    = @printf io "x -> (a+b)/2 + (b-a)/2*c*x/sqrt(1 + c^2 * x^2)\n"
Base.show(io::IO, t::Exponential)   = @printf io "x -> b + (1/c) * log(x-a)\n"

abstract AbstractParameter{T<:FloatingPoint}
abstract Parameter{T,U<:Transform} <: AbstractParameter{T}

typealias ParameterVector{T} Vector{AbstractParameter{T}}
typealias NullablePrior      Nullable{ContinuousUnivariateDistribution}

# do we need value bounds?
immutable UnscaledParameter{T,U} <: Parameter{T,U}
	key::Symbol
    value::T                    # transformed parameter value
    valuebounds::Interval{T}
    transbounds::Interval{T}
    prior::NullablePrior
    fixed::Bool
    description::String
end

immutable ScaledParameter{T,U} <: Parameter{T,U}
	key::Symbol
    value::T                    # scaled, transformed parameter value
    unscaledvalue::T			# unscaled, transformed parameter value
    scaling::Function
    valuebounds::Interval{T}
    transbounds::Interval{T}
    prior::NullablePrior
    fixed::Bool
    description::String
end

hasprior(p::Parameter) = !isnull(p.prior)

typealias NullableOrPrior Union(NullablePrior, ContinuousUnivariateDistribution)

# method to construct parameters
function parameter{T,U<:Transform}(key::Symbol,
	                    		   value::T,
	                    		   valuebounds::Interval{T} = (value,value),
	                    		   transbounds::Interval{T} = (value,value),
	                    		   transform::U             = Untransformed(),
	                    		   prior::NullableOrPrior   = NullablePrior();
	                    		   description::String      = "This variable is missing a description. Zac will not be pleased!",
								   fixed::Bool              = true,
								   scaling::Function        = identity)

	# ensure that we have a Nullable{Distribution}, if not construct one
	prior = !isa(prior,NullablePrior) ? NullablePrior(prior) : prior

	if scaling == identity
		return UnscaledParameter{T,U}(key, value, valuebounds, transbounds, prior, fixed, description)
	else
		return ScaledParameter{T,U}(key, scaling(value), value, scaling, valuebounds, transbounds, prior, fixed, description)
	end
end

# generate a parameter given a new value
function parameter{T,U}(p::UnscaledParameter{T,U}, newvalue::T)
	p.fixed && return p
	a,b = p.transbounds
	@assert a <= newvalue <= b "New value is out of bounds"
	UnscaledParameter{T,U}(p.key, newvalue, p.valuebounds, p.transbounds, p.prior, p.fixed, p.description)
end
function parameter{T,U}(p::ScaledParameter{T,U}, newvalue::T)
	p.fixed && return p
	a,b = p.transbounds
	@assert a <= newvalue <= b "New value is out of bounds"
	ScaledParameter{T,U}(p.key, p.scaling(newvalue), newvalue, p.scaling, p.valuebounds, p.transbounds, p.prior, p.fixed, p.description)
end

function Base.show{T,U}(io::IO, p::Parameter{T,U})
	@printf io "%s\n" typeof(p)
	@printf io ":%s\n%s\n"      p.key p.description
    @printf io "____________________________\n"
	@printf io "real value:        %+6f\n" toreal(p)
	!isa(U(),Untransformed) && @printf io "transformed value: %+6f\n" p.value

	if hasprior(p)
        @printf io "prior distribution:\n\t%s\n" get(p.prior)
    else
        @printf io "prior distribution:\n\t%s\n" "no prior"
	end

	@printf io "transformation:\n\t%s" U()
    @printf io "parameter is %s\n" p.fixed ? "fixed" : "not fixed"
end

# does anyone know what c does here?

# Untransformed
tomodel{T}(p::Parameter{T,Untransformed}, x::T) = x
toreal{T}(p::Parameter{T,Untransformed}, x::T = p.value) = x

# SquareRoot
function tomodel{T}(p::Parameter{T,SquareRoot}, x::T)
	(a,b), c = p.transbounds, one(T)
	(a+b)/2 + (b-a)/2*c*x/sqrt(1 + c^2 * x^2)
end
function toreal{T}(p::Parameter{T,SquareRoot}, x::T = p.value)
	(a,b), c = p.transbounds, one(T)
	cx = 2 * (x - (a+b)/2)/(b-a); (1/c)*cx/sqrt(1 - cx^2)
end

# Exponential
function tomodel{T}(p::Parameter{T,Exponential}, x::T)
	(a,b),c = p.transbounds,one(T)
	a + exp(x-b)
end
function toreal{T}(p::Parameter{T,Exponential}, x::T = p.value)
	(a,b),c = p.transbounds,one(T)
	b + (1/c) * log(x-a)
end

tomodel{T}(pvec::ParameterVector{T}) = map(tomodel, pvec)
toreal{T}(pvec::ParameterVector{T}, values::Vector{T}) = map(toreal, pvec, values)

# define operators to work on parameters

Base.convert{T<:Real}(::Type{T}, p::AbstractParameter) = convert(T,p.value)
Base.promote_rule{T<:Real,U<:Real}(::Type{AbstractParameter{T}}, ::Type{U}) = promote_rule(U,T)

Base.(:^)(p::AbstractParameter, x::Integer) = (^)(p.value, x)

for op in (:(Base.(:+)),
           :(Base.(:-)),
           :(Base.(:*)),
           :(Base.(:/)),
           :(Base.(:^)))

    @eval ($op)(p::AbstractParameter, q::AbstractParameter) = ($op)(p.value, q.value)
    @eval ($op)(p::AbstractParameter, x::Number)            = ($op)(p.value, x)
    @eval ($op)(x::Number, p::AbstractParameter)            = ($op)(x, p.value)
end

for f in (:(Base.exp),
          :(Base.log),
          :(Base.(:-)),
          :(Base.(:<)),
          :(Base.(:>)),
          :(Base.(:<=)),
          :(Base.(:>=)))

    @eval ($f)(p::AbstractParameter) = ($f)(p.value)
end

# this function is optimised for speed
function update!{T}(pvec::ParameterVector{T}, newvalues::Vector{T})
	@assert length(newvalues) == length(pvec) "Length of input vector (=$(length(newvalues))) must match length of parameter vector (=$(length(pvec)))"
   	map!(parameter, pvec, pvec, newvalues)
end
# define the non-mutating version like this because we need the type stability of map!
update{T}(pvec::ParameterVector{T}, newvalues::Vector{T}) = update!(copy(pvec), newvalues)

Distributions.pdf(p::AbstractParameter) = exp(logpdf(p))
Distributions.logpdf{T,U}(p::UnscaledParameter{T,U}) = logpdf(get(p.prior),p.value)
Distributions.logpdf{T,U}(p::ScaledParameter{T,U})   = logpdf(get(p.prior),p.unscaledvalue)

# this function is optimised for speed
function Distributions.logpdf{T}(pvec::ParameterVector{T})
	x = zero(T)
	@inbounds for i = 1:length(pvec)
        if hasprior(pvec[i])
    		x += logpdf(pvec[i])
        end
	end
	x
end

# calculate logpdf at new values, without needing to allocate a temporary array with update
function Distributions.logpdf{T}(pvec::ParameterVector{T}, newvalues::Vector{T})
    @assert length(newvalues) == length(pvec) "Length of input vector (=$(length(newvalues))) must match length of parameter vector (=$(length(pvec)))"

    x = zero(T)
    @inbounds for i = 1:length(pvec)
        if hasprior(pvec[i])
            x += logpdf(parameter(pvec[i], newvalues[i]))
        end
    end
    x
end

Distributions.pdf{T}(pvec::ParameterVector{T}) = exp(logpdf(pvec))
Distributions.pdf{T}(pvec::ParameterVector{T}, newvalues::Vector{T}) = exp(logpdf(pvec, newvalues))
