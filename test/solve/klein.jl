using DSGE, Test, JLD2, ModelConstructors, FileIO, LinearAlgebra, BenchmarkTools
import DSGE: klein_transition_matrices, n_model_states, n_backward_looking_states
import ModelConstructors: @test_matrix_approx_eq

regen = false

path    = dirname(@__FILE__)
ref_dir = joinpath(path, "reference")
ref_file = joinpath(ref_dir, "klein.jld2")

# --- Model setup -----------------------------------------------------------
# Mirrors the het_dsge solve path: steady state, then solve with the normalized
# Jacobian (m.testing = false), preserving nx across the toggle.
m = HetDSGE()
steadystate!(m)

nx_save = get_setting(m, :nx)
m.testing = false
m <= Setting(:nx, nx_save)

# --- Compute (always) ------------------------------------------------------
gx, hx, eu = klein(m)
# klein() returns (gx, hx, eu) = (TTT_jump, TTT_state, eu); callers then pass
# (TTT_state, TTT_jump) = (hx, gx) into klein_transition_matrices.
TTT, RRR = klein_transition_matrices(m, hx, gx)

# --- Regenerate or compare -------------------------------------------------
if regen
    isdir(ref_dir) || mkpath(ref_dir)
    JLD2.jldopen(ref_file, "w") do file
        write(file, "gx",  gx)
        write(file, "hx",  hx)
        write(file, "eu",  eu)
        write(file, "TTT", TTT)
        write(file, "RRR", RRR)
    end
    @info "Wrote klein reference outputs to $ref_file"
else
    saved_gx, saved_hx, saved_eu, saved_TTT, saved_RRR =
        JLD2.jldopen(ref_file, "r") do file
            read(file, "gx"), read(file, "hx"), read(file, "eu"),
            read(file, "TTT"), read(file, "RRR")
        end

    @testset "klein solution matches reference" begin
        # Success contract: callers (e.g. solve.jl, statespace_functions.jl) treat
        # a nonzero third return value as a solution failure.
        @test eu == 0
        @test eu == saved_eu
        @test @test_matrix_approx_eq saved_gx gx
        @test @test_matrix_approx_eq saved_hx hx
    end

    @testset "klein_transition_matrices matches reference" begin
        @test @test_matrix_approx_eq saved_TTT TTT
        @test @test_matrix_approx_eq saved_RRR RRR
    end
end

# --- Structural checks (no reference needed) -------------------------------
@testset "klein return shapes" begin
    NK = get_setting(m, :n_predetermined_variables)
    # gx maps the NK predetermined states to the non-predetermined block;
    # hx is the square state transition for the predetermined block.
    @test size(hx) == (NK, NK)
    @test size(gx, 2) == NK
    @test all(isfinite, gx)
    @test all(isfinite, hx)
end

@testset "klein_transition_matrices stacks the state space" begin
    nms  = n_model_states(m)
    nbls = n_backward_looking_states(m)

    @test size(TTT) == (nms, nms)
    @test size(RRR, 1) == nms

    # Backward-looking block carries the state transition hx.
    @test @test_matrix_approx_eq hx TTT[1:nbls, 1:nbls]
    # Jump block is loaded as gx*hx (time-t states -> time-t+1 jumps).
    @test @test_matrix_approx_eq gx * hx TTT[nbls+1:end, 1:nbls]
    # No model state maps forward from a jump: both right-hand blocks are zero.
    @test all(iszero, TTT[1:nbls, nbls+1:end])
    @test all(iszero, TTT[nbls+1:end, nbls+1:end])
end

# --- benchmarks ------------------------------------------------------
# Timings for the two solve routines. Set `run_benchmarks = false` to skip.
run_benchmarks = false
if run_benchmarks
    println("\nklein benchmarks (HetDSGE, nx = $(get_setting(m, :nx))):")
    print("  klein(m):                    "); @btime klein($m);
    print("  klein_transition_matrices:   "); @btime klein_transition_matrices($m, $hx, $gx);
end

nothing
