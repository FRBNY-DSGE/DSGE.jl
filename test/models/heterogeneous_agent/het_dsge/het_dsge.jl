using DSGE
using Test, BenchmarkTools
using JLD2

import DSGE: klein_transition_matrices, n_model_states, n_backward_looking_states

# What do you want to do?
check_steady_state = true

path = dirname(@__FILE__)

m = HetDSGE()

# Steady-state computation
if check_steady_state
    steadystate!(m)
    # @btime steadystate!(m)

    file = jldopen("$path/reference/steady_state.jld2", "r")
    saved_Rk   = read(file, "Rk")
    saved_ω    = read(file, "omega")
    saved_kl   = read(file, "kl")
    saved_k    = read(file, "k")
    saved_y    = read(file, "y")
    saved_T    = read(file, "T")

    saved_ell  = read(file, "ell")
    saved_c    = read(file, "c")
    saved_μ    = read(file, "mu")
    saved_β    = read(file, "beta")
    close(file)

    @testset "Check steady state outputs" begin
        @test saved_Rk ≈ m[:Rkstar].value
        @test saved_ω  ≈ m[:ωstar].value
        @test saved_kl ≈ m[:klstar].value
        @test saved_k  ≈ m[:kstar].value
        @test saved_y  ≈ m[:ystar].value
        @test saved_T  ≈ m[:Tstar].value

        @test saved_ell ≈ m[:lstar].value
        @test saved_c   ≈ m[:cstar].value
        @test saved_μ   ≈ m[:μstar].value
        # Tolerance of convergence of β is 1e-5
        @test abs(saved_β - m[:βstar].value) <= 1e-5
    end
end
