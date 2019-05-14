using DSGE
using Test, BenchmarkTools
using JLD2

# What do you want to do?
check_steady_state = true

path = dirname(@__FILE__)

het = HetDSGEGovDebt(testing_gamma = true)
het <= Setting(:steady_state_only, true)
steadystate!(het)

m   = RepDSGEGovDebt(het)

# Steady-state computation
if check_steady_state
    steadystate!(m, het)

    file = jldopen("$path/../../heterogeneous_agent/het_dsge_gov_debt/reference/steady_state.jld2", "r")
    saved_Rk   = read(file, "Rk")
    saved_ω    = read(file, "omega")
    saved_kl   = read(file, "kl")
    saved_k    = read(file, "k")
    saved_y    = read(file, "y")
    saved_T    = read(file, "T")

    close(file)

    file = jldopen("$path/reference/steady_state.jld2", "r")

    saved_c    = read(file, "c")
    saved_ell  = read(file, "l")
    saved_β    = read(file, "beta")

    close(file)

    @testset "Check steady state outputs" begin
        @test saved_Rk ≈ m[:Rkstar].value
        @test saved_ω  ≈ m[:ωstar].value
        @test saved_kl ≈ m[:klstar].value
        @test saved_k  ≈ m[:kstar].value
        @test saved_y  ≈ m[:ystar].value
        @test saved_T  ≈ m[:Tstar].value

        @test saved_c   ≈ m[:cstar].value
        @test saved_ell ≈ m[:lstar].value
        @test saved_β   ≈ m[:βstar].value
    end
end
