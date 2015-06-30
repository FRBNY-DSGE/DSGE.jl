using Base.Test
using Distributions

using DSGE: DistributionsExt, AbstractModel, M990



# src/AbstractModel.jl
function test_all()
    test_param()
    test_parameters()
    test_modelinds()
    println("All tests in AbstractModel.jl passed")
end



# src/abstractmodel/param.jl
function test_param()
    # Inner Param constructor
    α = Param(0.1596, false, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999))
    @test α.scalefunction == identity
    @test α.scaledvalue == 0.1596
    @test α.description == ""

    # Inner Param constructor, fixed = true
    α_fixed = Param(0.1596, true, (1e-5, 0.999), Normal(0.30, 0.05), 1, (1e-5, 0.999))
    @test isa(α_fixed.prior, PointMass)
    @test α_fixed.transformtype == 0
    @test α_fixed.transformbounds == (0.1596, 0.1596)

    # One-argument Param constructor
    δ = Param(0.025)
    @test δ.fixed
    @test δ.bounds == (0.025, 0.025)
    @test isa(δ.prior, PointMass)

    # Invalid transformtype
    @test_throws ErrorException α_bad = Param(0.1596, false, (1e-5, 0.999), Normal(0.30, 0.05), -1, (1e-5, 0.999))

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

    println("param.jl tests passed")
end



# src/abstractmodel/parameters.jl
function test_parameters()
    # TODO: test logprior

    println("parameters.jl tests passed")
end



# src/abstractmodel/modelinds.jl
function test_modelinds()
    # makedict
    dict = makedict(["alice", "bob"])
    @test dict["alice"] == 1
    @test dict["bob"] == 2
    @test_throws KeyError dict["carol"]

    println("modelinds.jl tests passed")
end



# src/abstractmodel/model.jl
# Model is an abstract type and there are no functions to test in this file
