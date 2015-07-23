using DSGE: DistributionsExt
using Distributions

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

# println("param.jl tests passed\n")

# Test Parameters type

# Parameters iterator
Θ = Parameters990(model_specifications(Model990))
lastparam = Param(0.0)
for α = Θ
    lastparam = α
end
@test isa(lastparam, Param)
@test lastparam.value == 0.0181

# logprior
priordensity = exp(logprior(Θ))
@test priordensity >= 0
@test priordensity <= 1

# println("parameters.jl tests passed\n")

# Test modelinds

# makedict
dict = makedict(["alice", "bob"])
@test dict["alice"] == 1
@test dict["bob"] == 2
@test_throws KeyError dict["carol"]

# makedict with start field
dict = makedict(["alice", "bob"]; start=4)
@test dict["alice"] == 5
@test dict["bob"] == 6

# println("modelinds.jl tests passed\n")
