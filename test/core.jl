using Distributions
import DSGE: Param, PointMass

# Test Param type

# Inner Param constructor
α = Param(0.1596, false, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999))
@test α.scalefunction == identity
@test α.scaledvalue == 0.1596
@test α.description == ""

# Inner Param constructor, fixed = true
α_fixed = Param(0.1596, true, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999))
@test isa(α_fixed.priordist, PointMass)
@test α_fixed.transformtype == 0

# One-argument Param constructor
δ = Param(0.025)
@test δ.fixed
@test δ.bounds == (0.025, 0.025)
@test isa(δ.priordist, PointMass)

# Invalid transformtype
@test_throws ErrorException α_bad = Param(0.1596, false, (1e-5, 0.999),
                                          Normal(0.30, 0.05), -1, (1e-5, 0.999))

# update! value
update!(α, 0.0)
@test α.value == 0.0
@test α.scaledvalue == 0.0
update!(δ, 0.0)
@test δ.value == 0.025
@test δ.scaledvalue == 0.025
update!(α, 0.1596)

# Arithmetic with parameters
@test convert(Float16, δ) == 0.025
@test promote_type(Param, Float16) == Float64
@test promote_type(Param, Int8) == Float64
@test δ + δ == 0.05
@test δ^2 == 0.025^2
@test -δ == -0.025
@test log(δ) == log(0.025)

# toreal and tomodel
cx = 2 * (α - 1/2)
@test_approx_eq_eps(toreal(α), cx / sqrt(1 - cx^2), 0.001)
@test toreal(δ) == 0.025


model = Model990()
lastparam = Param(0.0)
for α in model.par
    isa(α,Param) && (lastparam = α)
end
@test isa(lastparam, Param)
@test lastparam.value == 0.0181

# prior
priordensity = exp(prior(model))
@test priordensity >= 0
@test priordensity <= 1
