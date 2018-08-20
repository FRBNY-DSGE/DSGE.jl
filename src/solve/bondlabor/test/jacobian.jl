using Base.Test
using JLD

path = dirname(@__FILE__)

file = jldopen("$path/reference/jacobian.jld", "r")
saved_JJ  = read(file, "JJ")
close(file)

nx = get_setting(m, :nx)
ns = get_setting(m, :ns)

# we will always order things XP YP X Y
# convention is that capital letters generally refer to indices
MUP  = 1:nx*ns
ZP   = nx*ns+1
ELLP = nx*ns+2:2*nx*ns+1
RP   = 2*nx*ns+2

MU   = 2*nx*ns+3:3*nx*ns+2
Z    = 3*nx*ns+3
ELL  = 3*nx*ns+4:4*nx*ns+3
R    = 4*nx*ns+4

# create objects needed for solve.jl
F1 = 1:nx*ns # euler eqn
F2 = nx*ns+1:2*nx*ns # KF
F3 = 2*nx*ns+1:2*nx*ns+1 # mkt ckr
F4 = 2*nx*ns+2:2*nx*ns+2 # z

@testset "Check jacobian outputs" begin
    @testset "Euler Equation" begin
        @test saved_JJ[F1, ZP] ≈ JJ[F1, ZP]
        @test saved_JJ[F1, ELLP] ≈ JJ[F1, ELLP]
        @test saved_JJ[F1, RP] ≈ JJ[F1, RP]
        @test saved_JJ[F1, Z] ≈ JJ[F1, Z]
        @test saved_JJ[F1, ELL] ≈ JJ[F1, ELL]
        @test saved_JJ[F1, R] ≈ JJ[F1, R]
    end

    @testset "Kolmogorov Forward Equation" begin
        @test saved_JJ[F2, MUP] ≈ JJ[F2, MUP]
        @test saved_JJ[F2, MU] ≈ JJ[F2, MU]
        @test saved_JJ[F2, Z] ≈ JJ[F2, Z]
        @test saved_JJ[F2, ELL] ≈ JJ[F2, ELL]
        @test saved_JJ[F2, R] ≈ JJ[F2, R]
    end

    @testset "Market clearing condition" begin
        @test saved_JJ[F3, MU] ≈ JJ[F3, MU]
        @test saved_JJ[F3, Z] ≈ JJ[F3, Z]
        @test saved_JJ[F3, ELL] ≈ JJ[F3, ELL]
        @test saved_JJ[F3, R] ≈ JJ[F3, R]
    end

    @testset "Technology process" begin
        @test saved_JJ[F4, ZP] ≈ JJ[F4, ZP]
        @test saved_JJ[F4, Z] ≈ JJ[F4, Z]
    end
end
