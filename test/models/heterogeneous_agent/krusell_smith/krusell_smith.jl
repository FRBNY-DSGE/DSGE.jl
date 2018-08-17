using DSGE
using JLD
using Base.Test
import DSGE: jacobian

path = dirname(@__FILE__)

### Model
m = KrusellSmith()

### Steady State
file = jldopen("$path/reference/steady_state.jld", "r")
saved_l = read(file, "saved_l")
saved_c = read(file, "saved_c")
saved_μ = read(file, "saved_mu")
saved_K = read(file, "saved_K")
saved_KF = read(file, "saved_KF")
close(file)

@testset "Steady State" begin
    @test saved_l ≈ m[:lstar].value
    @test saved_c ≈ m[:cstar].value
    @test saved_μ ≈ m[:μstar].value
    @test saved_K ≈ m[:Kstar].value
    @test @test_matrix_approx_eq saved_KF m[:KFstar].value
end

### Jacobian
JJ       = jacobian(m)
saved_JJ = load("$path/reference/jacobian.jld", "saved_JJ")
nw       = load("$path/reference/jacobian.jld", "nw")

# we will always order things XP YP X Y
# convention is that capital letters generally refer to indices
# XP
LMP   = 1     :nw
LELLP = nw+1  :2*nw
KKP   = 2*nw+1

# YP
ZP    = 2*nw+2
MP    = 2*nw+3:3*nw+2
ELLP  = 3*nw+3:4*nw+2

# X
LM    = 4*nw+3:5*nw+2
LELL  = 5*nw+3:6*nw+2
KK    = 6*nw+3

# Y
Z     = 6*nw+4
M     = 6*nw+5:7*nw+4
ELL   = 7*nw+5:8*nw+4

# create objects needed for solve.jl
F1 = 1     :nw
F2 = nw+1  :2*nw
F3 = 2*nw+1:3*nw
F4 = 3*nw+1:4*nw
F5 = 4*nw+1:4*nw+1
F6 = 4*nw+2:4*nw+2

@testset "Jacobian" begin
    @testset "Euler Equation" begin
        @test @test_matrix_approx_eq saved_JJ[F1, KKP]  JJ[F1, KKP]
        @test @test_matrix_approx_eq saved_JJ[F1, ZP]   JJ[F1, ZP]
        @test @test_matrix_approx_eq saved_JJ[F1, ELLP] JJ[F1, ELLP]
        @test @test_matrix_approx_eq saved_JJ[F1, ELL]  JJ[F1, ELL]
    end

    @testset "Kolmogorov Forward Equation" begin
        @test @test_matrix_approx_eq saved_JJ[F2, LM]   JJ[F2, LM]
        @test @test_matrix_approx_eq saved_JJ[F2, LELL] JJ[F2, LELL]
        @test @test_matrix_approx_eq saved_JJ[F2, KK]   JJ[F2, KK]
        @test @test_matrix_approx_eq saved_JJ[F2, Z]    JJ[F2, Z]
        @test @test_matrix_approx_eq saved_JJ[F2, M]    JJ[F2, M]
    end

    @testset "LM(t+1) = M(t)" begin
        @test @test_matrix_approx_eq saved_JJ[F3, LMP]   JJ[F3, LMP]
        @test @test_matrix_approx_eq saved_JJ[F3, M]     JJ[F3, M]
    end

    @testset "LELL(t+1) = ELL(t)" begin
        @test @test_matrix_approx_eq saved_JJ[F4, LELLP]   JJ[F4, LELLP]
        @test @test_matrix_approx_eq saved_JJ[F4, ELL]     JJ[F4, ELL]
    end

    @testset "LOM K" begin
        @test @test_matrix_approx_eq saved_JJ[F5, ELL]   JJ[F5, ELL]
        @test @test_matrix_approx_eq saved_JJ[F5, M]     JJ[F5, M]
        @test @test_matrix_approx_eq saved_JJ[F5, KKP]   JJ[F5, KKP]
    end

    @testset "TFP" begin
        @test @test_matrix_approx_eq saved_JJ[F6, ZP]    JJ[F6, ZP]
        @test @test_matrix_approx_eq saved_JJ[F6, Z]     JJ[F6, Z]
    end
end

### Klein
TTT_jump, TTT_state = klein(m)

saved_gx = load("$path/reference/klein.jld", "gx")
saved_hx = load("$path/reference/klein.jld", "hx")

@testset "Solve" begin
    @test @test_matrix_approx_eq saved_gx TTT_jump
    @test @test_matrix_approx_eq saved_hx TTT_state
end
