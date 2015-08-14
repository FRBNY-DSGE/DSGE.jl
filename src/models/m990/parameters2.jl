using Distributions

type Param <: Number
    value::Float64
    scalefunction::Function
    scaledvalue::Float64
    fixed::Bool
    bounds::(Float64, Float64)
    priordist::Distribution
    transformtype::Int64
    transformbounds::(Float64, Float64)
    description::String
    texLabel::String
    
    function Param(value::Float64, fixed::Bool, bounds::(Float64, Float64),
                   priordist::Distribution, transformtype::Int64,
                   transformbounds::(Float64, Float64); scalefunction::Function = identity,
                   description::String = "",texLabel::String = "")
        if fixed
            priordist = PointMass(value)
            transformtype = 0
        end
        if transformtype != 0 && transformtype != 1 && transformtype != 2
            error("transformtype must be 0, 1, or 2")
        end
        (a, b) = transformbounds
        return new(value, scalefunction, scalefunction(value), fixed, bounds, priordist,
                   transformtype, transformbounds, description, texLabel)
    end
end

type Parameterisation{T <: AbstractModel}
	lookup::Dict{Symbol,Int}
	params::Vector
	description::String
end

function Parameterisation(T::Type{Model990})
	param_list = [:alp, :zeta_p, :iota_p, :del, :ups, :Bigphi, :s2, :h, :ppsi, :nu_l, :zeta_w, :iota_w, :law]
	lookup = Dict([(param_list[i],i) for i in 1:length(param_list)])

	params[lookup[:alp   ]] = Param(0.1596, false, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999))
	params[lookup[:zeta_p]] = Param(0.8940, false, (1e-5, 0.999), Beta(0.5, 0.1), 1, (1e-5, 0.999))
	params[lookup[:iota_p]] = Param(0.1865, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999))
	params[lookup[:del   ]] = Param(0.025)
	params[lookup[:ups   ]] = Param(1.000, true, (0., 10.), Gamma(1., 0.5), 2, (1e-5, 0.))
	params[lookup[:Bigphi]] = Param(1.1066, false, (1., 10.), Normal(1.25, 0.12), 2, (1.00, 10.00))
	params[lookup[:s2    ]] = Param(2.7314, false, (-15., 15.), Normal(4., 1.5), 0, (-15., 15.))
	params[lookup[:h     ]] = Param(0.5347, false, (1e-5, 0.999), Beta(0.7, 0.1), 1, (1e-5, 0.999))
	params[lookup[:ppsi  ]] = Param(0.6862, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999))
	params[lookup[:nu_l  ]] = Param(2.5975, false, (1e-5, 10.), Normal(2, 0.75), 2, (1e-5, 10.))
	params[lookup[:zeta_w]] = Param(0.9291, false, (1e-5, 0.999), Beta(0.5, 0.1), 1, (1e-5, 0.999))
	params[lookup[:iota_w]] = Param(0.2992, false, (1e-5, 0.999), Beta(0.5, 0.15), 1, (1e-5, 0.999))
	params[lookup[:law   ]] = Param(1.5)

	description = "this is some fancy stuff"
	Parameterisation{T}(dict, params, description)
end

Base.getindex(p::Parameterisation,ind) = p.params[p.lookup[ind]]
Base.keys(p::Parameterisation,ind) = keys(p.lookup)
