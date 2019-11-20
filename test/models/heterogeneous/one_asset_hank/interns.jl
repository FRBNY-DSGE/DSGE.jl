using DSGE, JLD2, FileIO, Distributions, Nullables, SparseArrays
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
    @test size(Γ0) == size(out["output"][1])
    @test size(Γ1) == size(out["output"][2])
    @test typeof(Γ0) == typeof(out["output"][1])
    @test typeof(Γ1) == typeof(out["output"][2])
    @test Γ0 == out["output"][1]
    @test Γ1 ≈ out["output"][2] # Note: this does not hold with perfect equality
    @test Ψ  == out["output"][3]
    @test Π  == out["output"][4]
    @test C  == out["output"][5]
end

# Test ouput of solve has not changed
#out = load("reference/solve_output_interns.jld2")
out = load("reference/one_asset_after_upgrade.jld2")
TTT, RRR, CCC, inverse_basis = solve(m)
@testset "solve(m) output" begin
    @test TTT           = out["TTT"]
    @test RRR           = out["RRR"]
    @test CCC           = out["CCC"]
    @test inverse_basis = out["inverse_basis"]
end
