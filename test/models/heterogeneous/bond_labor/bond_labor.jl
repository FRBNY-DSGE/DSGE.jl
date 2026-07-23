using DSGE
using Test, BenchmarkTools
using JLD2

path = dirname(@__FILE__)

run_benchmarks = false

m = BondLabor()

# Steady-state computation
steadystate!(m)
run_benchmarks && @btime steadystate!(m)

file = JLD2.jldopen("$path/reference/steady_state.jld2", "r")
saved_ell  = read(file, "ell")
saved_c    = read(file, "c")
saved_η    = read(file, "eta")
saved_μ    = read(file, "mu")
saved_β    = read(file, "beta")
saved_χss  = read(file, "chi_ss")
close(file)

@testset "Check steady state outputs" begin
    @test saved_ell ≈ m[:lstar].value
    @test saved_c   ≈ m[:cstar].value
    @test saved_η   ≈ m[:ηstar].value
    @test saved_μ   ≈ m[:μstar].value
    @test saved_β   ≈ m[:βstar].value
    @test saved_χss ≈ m[:χstar].value
end

# Jacobian computation
m.testing = true        # So that it will test against the unnormalized Jacobian
JJ = DSGE.jacobian(m)
run_benchmarks && @btime JJ = DSGE.jacobian(m)

file = JLD2.jldopen("$path/reference/jacobian.jld2", "r")
saved_JJ  = read(file, "JJ")
close(file)

nx = DSGE.get_setting(m, :nx)
ns = DSGE.get_setting(m, :ns)

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
run_benchmarks && @btime klein(m)

@JLD2.load "$path/reference/solve.jld2" saved_gx saved_hx

@testset "Check solve outputs" begin
    @test saved_gx  ≈ gx
    @test saved_hx  ≈ hx
end

# State-space transition matrices (klein returns TTT_jump, TTT_state = gx, hx)
TTT, RRR = DSGE.klein_transition_matrices(m, hx, gx)
run_benchmarks && @btime DSGE.klein_transition_matrices(m, hx, gx)

CCC = zeros(DSGE.n_model_states(m))

# Shock loading
RRR_shock = DSGE.shock_loading(m, gx)
run_benchmarks && @btime DSGE.shock_loading(m, gx)

@testset "Check shock loading" begin
    nb   = DSGE.n_backward_looking_states(m)
    exo  = m.exogenous_shocks
    endo = m.endogenous_states
    # Technology shock loads unity onto the z′ state...
    @test RRR_shock[endo[:z′_t], exo[:z_sh]] ≈ ones(1)
    # ...and the jump block is the jump policy applied to the state loading
    @test RRR_shock[nb+1:end, :] ≈ gx * RRR_shock[1:nb, :]
    # klein_transition_matrices must reuse the same shock loading
    @test RRR ≈ RRR_shock
end

# Measurement equation
meas = DSGE.measurement(m, TTT, gx, RRR, CCC)
run_benchmarks && @btime DSGE.measurement(m, TTT, gx, RRR, CCC)

@testset "Check measurement equation" begin
    obs = m.observables
    exo = m.exogenous_shocks
    nb  = DSGE.n_backward_looking_states(m)
    # Dimensions
    @test size(meas.ZZ) == (DSGE.n_observables(m), DSGE.n_model_states(m))
    @test length(meas.DD) == DSGE.n_observables(m)
    # Innovation variance and measurement error are set by construction
    @test meas.QQ[exo[:z_sh], exo[:z_sh]] ≈ m[:σ_z].value^2
    @test meas.EE[obs[:obs_gdp], obs[:obs_gdp]] ≈ m[:e_y].value
    # CCC = 0 ⟹ no intercept adjustment
    @test all(meas.DD .== 0)
    # GDP loads only on states, so the jump columns of ZZ vanish
    @test all(meas.ZZ[:, nb+1:end] .== 0)
end
