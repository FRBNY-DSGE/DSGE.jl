using Compat, Distributions, Debug
import DSGE: PointMass

typealias Interval{T} @compat Tuple{T,T}

# define all the kinds of transformations we make
abstract Transform
type Untransformed <: Transform end
type SquareRoot    <: Transform end
type Exponential   <: Transform end

Base.show(io::IO, t::Untransformed) = @printf io "x -> x\n"
Base.show(io::IO, t::SquareRoot)    = @printf io "x -> (a+b)/2 + (b-a)/2*c*x/sqrt(1 + c^2 * x^2)\n"
Base.show(io::IO, t::Exponential)   = @printf io "x -> b + (1/c) * log(x-a)\n"

abstract AbstractParameter{T<:Real} # made this real for conversions
abstract Parameter{T,U<:Transform} <: AbstractParameter{T}

typealias ParameterVector{T} Vector{AbstractParameter{T}}
typealias NullablePrior      Nullable{ContinuousUnivariateDistribution}

type UnscaledParameter{T,U} <: Parameter{T,U}
    key::Symbol
    value::T                             # parameter value in model space
    valuebounds::Interval{T}             # bounds of parameter value
    transform_parameters::Interval{T}    # parameters for transformation
    transform::U                # transformation we use to go between model space and real line for csminwel
    prior::NullablePrior        # prior distribution
    fixed::Bool                 # is this parameter fixed at some value, or do we estimate it?
    description::AbstractString 
    texLabel::AbstractString    # LaTeX label for printing
end

type ScaledParameter{T,U} <: Parameter{T,U}
    key::Symbol
    value::T                    # unscaled parameter value in model space
    scaledvalue::T		# scaled parameter value in model space
    valuebounds::Interval{T}    # 
    transform_parameters::Interval{T}    # 
    transform::U                # we only use transformed values for csminwel.
                                # tomodel/toreal takes care of the conversion.
    prior::NullablePrior
    fixed::Bool
    scaling::Function
    description::AbstractString
    texLabel::AbstractString
end

"""
type SteadyStateParameter{T} <: AbstractParameter{T}

Defines a type for the model's steady-state parameters, stored in `m.steady_state`. Their values are calculated and set by `steadystate!(m)`, rather than being estimated directly.  `SteadyStateParameter`s have the following fields:

-`key::Symbol`: Key for referencing this parameter

-`value::T`: The parameter's steady-state value

-`description::AbstractString`: Short description of the parameter's economic significance

-`texLabel::AbstractString`: for LaTeX printing

Unlike true model parameters (for which we must find modal values using csminwel), these values need not be transformed from the model space to the real line. Thus, they do not have `transform` fields. 
"""
type SteadyStateParameter{T} <: AbstractParameter{T}
    key::Symbol
    value::T                    
    description::AbstractString
    texLabel::AbstractString
end

hasprior(p::Parameter) = !isnull(p.prior)

typealias NullableOrPrior Union(NullablePrior, ContinuousUnivariateDistribution)

# We want to use values from UnscaledParameters and
# SteadyStateParameters in computation, so we alias their union here.
typealias UnscaledOrSteadyState Union(UnscaledParameter, SteadyStateParameter)

# method to construct parameters
function parameter{T,U<:Transform}(key::Symbol,
                                   value::T,
                                   valuebounds::Interval{T} = (value,value),
                                   transform_parameters::Interval{T} = (value,value),
                                   transform::U             = Untransformed(),
                                   prior::NullableOrPrior   = NullablePrior();
                                   fixed::Bool              = true,
                                   scaling::Function        = identity,
                                   description::AbstractString      = "This variable is missing a description.",
                                   texLabel::AbstractString = "")

    
    # If fixed=true, force bounds to match and prior to be
    # PointMass. We need to define new variable names here because of lexical
    # scoping.

    ret_valuebounds = valuebounds
    ret_transform_parameters = transform_parameters
    ret_prior = prior

    if fixed
        ret_transform_parameters = (value,value)  # value is transformed already       
        ret_prior = PointMass(value)

        if isa(transform, Untransformed)
            ret_valuebounds = (value,value)
        end
    else
        ret_transform_parameters = transform_parameters
    end
    
    # ensure that we have a Nullable{Distribution}, if not construct one
    ret_prior = !isa(ret_prior,NullablePrior) ? NullablePrior(ret_prior) : ret_prior

    if scaling == identity
        return UnscaledParameter{T,U}(key, value, ret_valuebounds, ret_transform_parameters, transform,
                                      ret_prior, fixed, description, texLabel)
    else
        return ScaledParameter{T,U}(key, value, scaling(value), ret_valuebounds, ret_transform_parameters, transform,
                                    ret_prior, fixed, scaling, description, texLabel)
    end
end

# construct steady-state parameters
function SteadyStateParameter{T<:Real}(key::Symbol,
                                       value::T;
                                       description::AbstractString = "No description provided.",
                                       texLabel::AbstractString = "")

    return SteadyStateParameter(key, value, description, texLabel)
end


# generate a parameter given a new value 
function parameter{T<:Real,U<:Transform}(p::UnscaledParameter{T,U}, newvalue::T)
    p.fixed && return p  # if the parameter is fixed, don't change its value
    a,b = p.valuebounds  
    @assert a <= newvalue <= b "New value is out of bounds"
    UnscaledParameter{T,U}(p.key, newvalue, p.valuebounds, p.transform_parameters, p.prior, p.fixed, p.description)
end
function parameter{T<:Real,U<:Transform}(p::ScaledParameter{T,U}, newvalue::T)
    p.fixed && return p
    a,b = p.valuebounds  
    @assert a <= newvalue <= b "New value is out of bounds"
    ScaledParameter{T,U}(p.key, newvalue, p.scaling(newvalue), p.scaling, p.valuebounds, p.transform_parameters, p.prior, p.fixed, p.description)
end

function Base.show{T,U}(io::IO, p::Parameter{T,U})
    @printf io "%s\n" typeof(p)
    @printf io "(:%s)\n%s\n"      p.key p.description
    @printf io "LaTeX label: %s\n"     p.texLabel
    @printf io "-----------------------------\n"
    #@printf io "real value:        %+6f\n" toreal(p)
    @printf io "unscaled, untransformed value:        %+6f\n" p.value
    isa(p,ScaledParameter) && @printf "scaled, untransformed value:        %+6f\n" p.scaledvalue
    #!isa(U(),Untransformed) && @printf io "transformed value: %+6f\n" p.value
    
    if hasprior(p)
        @printf io "prior distribution:\n\t%s\n" get(p.prior)
    else
        @printf io "prior distribution:\n\t%s\n" "no prior"
    end

    @printf io "transformation for csminwel:\n\t%s" U()
    @printf io "parameter is %s\n" p.fixed ? "fixed" : "not fixed"
end

function Base.show{T}(io::IO, p::SteadyStateParameter{T})
    @printf io "%s\n" typeof(p)
    @printf io "(:%s)\n%s\n"      p.key p.description
    @printf io "LaTeX label: %s\n"     p.texLabel
    @printf io "-----------------------------\n"
    @printf io "value:        %+6f\n" p.value
end


# does anyone know what c does here?

# Untransformed
tomodel{T}(p::Parameter{T,Untransformed}, x::T) = x
toreal{T}(p::Parameter{T,Untransformed}, x::T = p.value) = x

# SquareRoot
function tomodel{T}(p::Parameter{T,SquareRoot}, x::T)
    (a,b), c = p.transform_parameters, one(T)
    (a+b)/2 + (b-a)/2*c*x/sqrt(1 + c^2 * x^2)
end
function toreal{T}(p::Parameter{T,SquareRoot}, x::T = p.value)
    (a,b), c = p.transform_parameters, one(T)
    cx = 2 * (x - (a+b)/2)/(b-a)
    (1/c)*cx/sqrt(1 - cx^2)
end

# Exponential
function tomodel{T}(p::Parameter{T,Exponential}, x::T)
    (a,b),c = p.transform_parameters,one(T)
    a + exp(c*(x-b))
end
function toreal{T}(p::Parameter{T,Exponential}, x::T = p.value)
    (a,b),c = p.transform_parameters,one(T)
    b + (1/c) * log(x-a)
end

tomodel{T}(pvec::ParameterVector{T}) = map(tomodel, pvec)
toreal{T}(pvec::ParameterVector{T}, values::Vector{T}) = map(toreal, pvec, values)


# define operators to work on parameters

# TODO: do we also want to convert p to type AbstractParameter{T}? Seems so.
Base.convert{T<:Real}(::Type{T}, p::UnscaledParameter)  = convert(T,p.value)
Base.convert{T<:Real}(::Type{T}, p::ScaledParameter)    = convert(T,p.scaledvalue)  
Base.convert{T<:Real}(::Type{T}, p::SteadyStateParameter)  = convert(T,p.value)

Base.promote_rule{T<:Real,U<:Real}(::Type{AbstractParameter{T}}, ::Type{U}) = promote_rule(T,U)

for op in (:(Base.(:+)),
           :(Base.(:-)),
           :(Base.(:*)),
           :(Base.(:/)),
           :(Base.(:^)))

    @eval ($op)(p::UnscaledOrSteadyState, q::UnscaledOrSteadyState) = ($op)(p.value, q.value)
    @eval ($op)(p::UnscaledOrSteadyState, x::Number)            = ($op)(p.value, x)
    @eval ($op)(x::Number, p::UnscaledOrSteadyState)            = ($op)(x, p.value)

    @eval ($op)(p::ScaledParameter, q::ScaledParameter) = ($op)(p.scaledvalue, q.scaledvalue)
    @eval ($op)(p::ScaledParameter, x::Number)            = ($op)(p.scaledvalue, x)
    @eval ($op)(x::Number, p::ScaledParameter)            = ($op)(x, p.scaledvalue)

    @eval ($op)(p::ScaledParameter, q::UnscaledOrSteadyState) = ($op)(p.scaledvalue, q.value)
    @eval ($op)(p::UnscaledOrSteadyState, q::ScaledParameter) = ($op)(p.value, q.scaledvalue)
end

for f in (:(Base.exp),
          :(Base.log),
          :(Base.(:-)),
          :(Base.(:<)),
          :(Base.(:>)),
          :(Base.(:<=)),
          :(Base.(:>=)))

    @eval ($f)(p::UnscaledOrSteadyState) = ($f)(p.value)
    @eval ($f)(p::ScaledParameter) = ($f)(p.scaledvalue)
end

# this function is optimised for speed
function update!{T}(pvec::ParameterVector{T}, newvalues::Vector{T})
    @assert length(newvalues) == length(pvec) "Length of input vector (=$(length(newvalues))) must match length of parameter vector (=$(length(pvec)))"
   	map!(parameter, pvec, pvec, newvalues)
end
# define the non-mutating version like this because we need the type stability of map!
update{T}(pvec::ParameterVector{T}, newvalues::Vector{T}) = update!(copy(pvec), newvalues)

Distributions.pdf(p::AbstractParameter) = exp(logpdf(p))
Distributions.logpdf{T,U}(p::Parameter{T,U}) = logpdf(get(p.prior),p.value) # we want the unscaled value for ScaledParameters


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
