using DSGE, StateSpaceRoutines
using JLD
using Base.Test
import DSGE: jacobian, n_observables, n_shocks_exogenous, n_backward_looking_states, n_model_states

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

### Filtering

# Setup
N    = 200
m    = KrusellSmith()
endo = m.endogenous_states

# Model Solution/Transition Equation
TTT, RRR, CCC = solve(m)

# Measurement Equation
meas = measurement(m, TTT, RRR, CCC)
ZZ  = meas[:ZZ]
DD  = meas[:DD]
EE  = fill(0.1, (1,1))
QQ  = meas[:QQ]

# Generate measurement errors and shocks
srand(42)
u_t     = EE*randn(n_observables(m), N)
ε_t     = QQ*randn(n_shocks_exogenous(m), N)

####################################
# Simulate states forward N periods
####################################
s_t = zeros(n_model_states(m), N)
for t in 1:N
    # Initialize with shock around steady state
    if t == 1
        s_t[:, t] = RRR*ε_t[:, t]
    else
        s_t[:, t] = TTT*s_t[:, t-1] + RRR*ε_t[:, t]
    end
end

# Simulation in terms of grid values
simulated_log_gdp = vec(s_t[endo[:z′_t], :] + (m[:α]/m[:Kstar])*s_t[endo[:K′_t], :])
meas_log_gdp      = simulated_log_gdp + u_t[1, :]

##########
# Filter!
##########

# Measurement equation
data = ZZ*s_t + u_t   # Should be the same as meas_log_gdp

# Initialize with 0 initial condition
s_0     = zeros(n_model_states(m))
P_0     = RRR*QQ*RRR'

# Testing the simulation and filter setup
file = jldopen("$path/reference/simulate_and_filter.jld", "r")
saved_u_t  = read(file, "u_t")
saved_ε_t  = read(file, "eps_t")
saved_s_t  = read(file, "s_t")
saved_s_jump_t = read(file, "s_jump_t")
saved_simulated_log_gdp = read(file, "simulated_log_gdp")
saved_meas_log_gdp = read(file, "meas_log_gdp")
saved_data = read(file, "data")
saved_s_0  = read(file, "s_0")
saved_P_0  = read(file, "P_0")
close(file)

@testset "Simulate and Filter" begin
    @test saved_u_t ≈ u_t
    @test saved_ε_t ≈ ε_t
    @test saved_s_t ≈ s_t[1:n_backward_looking_states(m), :]
    @test saved_s_jump_t ≈ s_t[n_backward_looking_states(m)+1:end, :]
    @test saved_simulated_log_gdp ≈ simulated_log_gdp
    @test saved_meas_log_gdp ≈ meas_log_gdp
    @test saved_data ≈ data
    @test saved_s_0 ≈ s_0[1:n_backward_looking_states(m)]
    @test saved_P_0 ≈ P_0[1:n_backward_looking_states(m), 1:n_backward_looking_states(m)]
end

# Testing the filter inputs against the outputs
EE_zero = zeros(1, 1)
reference_output = load("$path/reference/filter_inputs_output.jld")
test = kalman_filter(reference_output["y"], TTT, RRR, CCC, QQ, ZZ,
                     DD, EE_zero, s_0, P_0)

@testset "Check Filter Output" begin
    @test test[1] ≈ reference_output["loglik"]
end
