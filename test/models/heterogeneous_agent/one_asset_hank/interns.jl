using DSGE, JLD2, FileIO, Distributions, Nullables
import Test: @test, @testset, @test
import DataStructures: OrderedDict

### Model
m = OneAssetHANK()
# Test ouput of steadystate! has not changed
#=
out = load("reference/new_steadystate_output.jld2")
@testset "steadystate!(m) output" begin
    steadystate!(m)
    N = length(m.parameters)
    if (m.steady_state[N+1].key == :V_ss)
        N += 1
    end
    for param in out["parameters"]
        @test m.parameters[m.keys[param.key]].value == param.value
    end
    for ss_param in out["steady_state"]
        @test m.steady_state[m.keys[ss_param.key] - N].value == ss_param.value
    end
    @test m.keys == out["keys"]
end
=#

# Test ouput of eqcond has not changed
out = load("reference/eqcond_output.jld2")
Γ0, Γ1, Ψ, Π, C = eqcond(m)
@testset "eqcond(m) output" begin
    @test size(Γ0) == size(out["Γ0"])
    @test size(Γ1) == size(out["Γ1"])
    @test typeof(Γ0) == typeof(out["Γ0"])
    @test typeof(Γ1) == typeof(out["Γ1"])
    @test Γ0 == out["Γ0"]
    @test Γ1 == out["Γ1"]
    @test Ψ  == out["Ψ"]
    @test Π  == out["Π"]
    @test C  == out["C"]
end

# Test ouput of solve has not changed
out = load("reference/new_solve_output.jld2")
TTT, RRR, CCC, inverse_basis = solve(m)
@testset "solve(m) output" begin
    @test TTT           ≈ out["TTT"]
    @test RRR           ≈ out["RRR"]
    @test CCC           ≈ out["CCC"]
    @test inverse_basis ≈ out["inverse_basis"]
end
