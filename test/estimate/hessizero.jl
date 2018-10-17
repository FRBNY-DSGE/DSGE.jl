using DSGE
using Test

# Test `hessizero` in context of Rosenbrock function
function rosenbrock(x::Vector)
    a = 100
    b = 1
    return a*(x[2]-x[1]^2.0)^2.0 + b*(1-x[1])^2.0
end

function rosenbrock_hessian(x::Vector)
    H = zeros(eltype(x), 2, 2)
    H[1,1] = 2.0 - 400.0 * x[2] + 1200.0 * x[1]^2.0
    H[1,2] = -400.0 * x[1]
    H[2,1] = -400.0 * x[1]
    H[2,2] = 200.0
    return H
end

# At the min, ensure no negatives in diagonal
x0 = [1.0, 1.0]
hessian_expected = rosenbrock_hessian(x0)
hessian, = DSGE.hessizero(rosenbrock, x0; check_neg_diag=true)
@testset "Check valid hessian at the minimum" begin
    @test hessian_expected ≈ hessian atol=0.01
    #@test_matrix_approx_eq hessian_expected hessian
end

# Not at the min (indeed, we are at max), ensure throws error
x1 = [1.0, 1.0]
rosenbrock_neg(x) = -rosenbrock(x)
@testset "Throw error for hessian calculation not at the min" begin
    @test_throws Exception DSGE.hessizero(rosenbrock_neg, x1; check_neg_diag=true)
end

# Not at the min, check matches closed form
x1 = Vector[[0.5, 1.5], [-1.0, -1.0]]
@testset "Check closed-form of the hessian" begin
    for x in x1
        local hessian_expected = rosenbrock_hessian(x)
        local hessian, = DSGE.hessizero(rosenbrock, x)
        @test @test_matrix_approx_eq hessian_expected hessian
    end
end
