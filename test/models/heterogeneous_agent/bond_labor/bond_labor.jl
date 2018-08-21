using DSGE
using Base.Test, BenchmarkTools
using JLD

path = dirname(@__FILE__)

m = BondLabor()

# Steady-state computation
steadystate!(m)
@btime steadystate!(m)

file = jldopen("$path/reference/steady_state.jld", "r")
saved_ell  = read(file, "ell")
saved_c    = read(file, "c")
saved_η    = read(file, "eta")
saved_μ    = read(file, "mu")
saved_β    = read(file, "beta")
close(file)

@testset "Check steady state outputs" begin
    @test saved_ell ≈ m[:lstar].value
    @test saved_c   ≈ m[:cstar].value
    @test saved_η   ≈ m[:ηstar].value
    @test saved_μ   ≈ m[:μstar].value
    @test saved_β   ≈ m[:βstar].value
end

# Jacobian computation
m.testing = true        # So that it will test against the unnormalized Jacobian
JJ = DSGE.jacobian(m)
@btime JJ = DSGE.jacobian(m)

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

# Solve
m.testing = false      # So the Jacobian will be normalized within the klein solution
gx, hx = klein(m)
@btime klein(m)

file = jldopen("$path/reference/solve.jld", "r")
# saved_gx2  = read(file, "gx2")
# saved_hx2  = read(file, "hx2")
saved_gx   = read(file, "gx")
saved_hx   = read(file, "hx")
close(file)

@testset "Check solve outputs" begin
    # @test saved_gx2 ≈ gx2
    # @test saved_hx2 ≈ hx2
    @test saved_gx  ≈ gx
    @test saved_hx  ≈ hx
end
