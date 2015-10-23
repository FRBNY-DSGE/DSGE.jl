using Distributions, Compat
import DSGE: PointMass

# Test Parameter type

# UnscaledParameter, fixed=false
α =  parameter(:α, 0.1596, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Normal(0.30, 0.05), fixed=false)
@test isa(α, UnscaledParameter)
@test α.key == :α
@test isa(α.prior.value, Normal)
@test α.prior.value.μ == 0.3
@test α.description == "This variable is missing a description."
@test α.texLabel == ""

# UnscaledParameter, fixed = true
α_fixed =  parameter(:α_fixed, 0.1596, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Normal(0.30, 0.05), fixed=true)
@test α_fixed.transbounds == (0.1596,0.1596)
@test isa(α_fixed.prior.value, PointMass) 
@test isa(α_fixed.transform, SquareRoot)


# Fixed UnscaledParameter, minimal constructor
δ = parameter(:δ, 0.025)
@test δ.fixed
@test δ.transbounds == (0.025, 0.025)
@test δ.valuebounds == (0.025, 0.025)
println(δ.prior.value)
@test isa(δ.prior.value, PointMass)

# Scaled parameter
β = parameter(:β, 0.1402, (1e-5, 10.), (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.25, 0.1), fixed=false,  scaling = x -> (1 + x/100)\1, description="β: Discount rate.", texLabel="\\beta ")
@test isa(β, ScaledParameter)
@test isa(β.prior.value, Gamma)
@test isa(β.transform, DSGE.Exponential)

# Invalid transform
@test_throws UndefVarError α_bad = parameter(:α, 0.1596, (1e-5, 0.999), (1e-5, 0.999),
                                              InvalidTransform(), Normal(0.30, 0.05), fixed=false)

## # update! value
## update!(α, 0.0)
## @test α.value == 0.0

## # update! fixed parameter shouldn't change its value
## update!(δ, 0.0)
## @test δ.value == 0.025

## # update! scaled parameter
## update!(β, 0.0)
## @test α.unscaledvalue == 0.0


## update!(α, 0.1596)

# Arithmetic with parameters
@test promote_type(AbstractParameter{Float64}, Float16) == Float64
@test promote_type(AbstractParameter{Float64}, Int8) == Float64
## @test promote_rule(AbstractParameter{Float64}, Float16) == Float64
## @test promote_rule(AbstractParameter{Float64}, Int8) == Float64
@test δ + δ == 0.05
@test δ^2 == 0.025^2
@test -δ == -0.025
@test log(δ) == log(0.025)

# toreal and tomodel
cx = 2 * (α - 1/2)
@test_approx_eq_eps(toreal(α), cx / sqrt(1 - cx^2), 0.001)
@test toreal(δ) == 0.025


model = Model990()
lastparam = parameter(:p, 0.0)
for α in model.parameters
    isa(α, Parameter) && (lastparam = α)
end
@test isa(lastparam, Parameter)
@test lastparam.value == 0.0181

# prior
priordensity = exp(prior(model))
@test priordensity >= 0
@test priordensity <= 1
