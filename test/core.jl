using Distributions, Compat


# Test Parameter type

# UnscaledParameter, fixed=false
α =  parameter(:α, 0.1596, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Normal(0.30, 0.05), fixed=false)
@test isa(α, UnscaledParameter)
@test α.key == :α
@test isa(α.prior.value, Normal)
@test α.prior.value.μ == 0.3
@test α.description == ""
@test α.texLabel == ""

# UnscaledParameter, fixed = true
α_fixed =  parameter(:α_fixed, 0.1596, (1e-5, 0.999), (1e-5, 0.999), SquareRoot(), Normal(0.30, 0.05), fixed=true)
@test α_fixed.transform_parameterization == (0.1596,0.1596)
@test isa(α_fixed.transform, SquareRoot)


# Fixed UnscaledParameter, minimal constructor
δ = parameter(:δ, 0.025)
@test δ.fixed
@test δ.transform_parameterization == (0.025, 0.025)
@test δ.valuebounds == (0.025, 0.025)

# Scaled parameter
β = parameter(:β, 0.1402, (1e-5, 10.), (1e-5, 10.), DSGE.Exponential(), GammaAlt(0.25, 0.1), fixed=false,  scaling = x -> (1 + x/100)\1, description="β: Discount rate.", texLabel="\\beta ")
@test isa(β, ScaledParameter)
@test isa(β.prior.value, Gamma)
@test isa(β.transform, DSGE.Exponential)

# Invalid transform
@test_throws UndefVarError α_bad = parameter(:α, 0.1596, (1e-5, 0.999), (1e-5, 0.999),
                                              InvalidTransform(), Normal(0.30, 0.05), fixed=false)


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

# settings
# settings - boolean, string, and number. adding to model. overwriting. filestrings. testing/not testing.
num_mh_blocks = Setting(:num_mh_blocks, 22) # short constructor
reoptimize = Setting(:reoptimize, false)
data_vintage = Setting(:data_vintage, "REF", true, "vint", "Date of data") # full constructor
@test promote_rule(Setting{Float64}, Float16) == Float64
@test promote_rule(Setting{Bool}, Bool) == Bool
@test promote_rule(Setting{ASCIIString}, AbstractString) == UTF8String
@test convert(Int64, num_mh_blocks) == 22
@test convert(ASCIIString, data_vintage) == "REF"


@test get_setting(model, :num_mh_blocks) == model.settings[:num_mh_blocks].value
model.testing = true
@test get_setting(model, :num_mh_blocks) == model.test_settings[:num_mh_blocks].value
@test modelstring(model) == "_test"

model.testing = false
model <= Setting(:num_mh_blocks, 5, true, "mhbk", "Number of blocks for Metropolis-Hastings")
@test model.settings[:num_mh_blocks].value == 5
@test ismatch(r"^\s*vint=(\d{6})_mhbk=5", modelstring(model))


nothing

