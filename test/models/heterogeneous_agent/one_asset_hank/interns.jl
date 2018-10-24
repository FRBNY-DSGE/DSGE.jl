using DSGE, JLD
import Base.Test: @test, @testset
import DataStructures: OrderedDict

### Model
m = OneAssetHANK()
# Test ouput of steadystate! has not changed
out = load("reference/steadystate_output_interns.jld")
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

# Test ouput of eqcond has not changed
out = load("reference/eqcond_output_interns.jld")
@testset "eqcond(m) output" begin
    Γ0, Γ1, Ψ, Π, C = eqcond(m)
    @test size(Γ0) == size(out["output"][1])
    @test size(Γ1) == size(out["output"][2])
    @test typeof(Γ0) == typeof(out["output"][1])
    @test typeof(Γ1) == typeof(out["output"][2])
    @test Γ0 == out["output"][1]
    @test Γ1 == out["output"][2]
    @test Ψ  == out["output"][3]
    @test Π  == out["output"][4]
    @test C  == out["output"][5]
end

# Test ouput of solve has not changed
out = load("reference/solve_output_interns.jld")
@testset "solve(m) output" begin
    TTT, RRR, CCC, inverse_basis = solve(m)
    @test TTT           == out["TTT"]
    @test RRR           == out["RRR"]
    @test CCC           == out["CCC"]
    @test inverse_basis == out["inverse_basis"]
end
