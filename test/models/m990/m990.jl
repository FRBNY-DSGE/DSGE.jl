using DSGE
using HDF5, Test, Distributions, JLD2

path = dirname(@__FILE__)

### Model
sub = "ss2"
model = Model990(sub)

### Parameters

@testset "Test parameter bounds checking" begin
    for θ in model.parameters
        !isa(θ,AbstractParameter) && error()
        θ.fixed && continue

        (left, right) = θ.valuebounds
        @test left < θ.value < right
    end
end

### Model indices

@testset "Checking model indices" begin
    # Endogenous states
    endo = model.endogenous_states
    @test length(endo) == 60
    @test endo[:Ez_t] == 60

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 16
    @test exo[:corepce_sh] == 16

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 13
    @test ex[:Erk_f_sh] == 13

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 60
    @test eq[:eq_Ez] == 60

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 12
    @test endo_new[:y_t1] == 61

    # Observables
    obs = model.observables
    @test length(obs) == 12
    @test obs[:obs_tfp] == 12
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Matrices are of expected dimensions
@testset "Check eqcond and state-space dimensions" begin
    @test size(Γ0) == (60, 60)
    @test size(Γ1) == (60, 60)
    @test size(C) == (60,)
    @test size(Ψ) == (60, 16)
    @test size(Π) == (60, 13)
end

# Check output matrices against reference output (ϵ = 1e-4)
h5 = h5open("$path/eqcond.h5")
Γ0_ref = read(h5, "G0")
Γ1_ref = read(h5, "G1")
C_ref  = reshape(read(h5, "C"), 60, 1)
Ψ_ref  = read(h5, "PSI")
Π_ref  = read(h5, "PIE")
close(h5)

@testset "Checking eqconds against reference output" begin
    @test @test_matrix_approx_eq Γ0_ref Γ0
    @test @test_matrix_approx_eq Γ1_ref Γ1
    @test @test_matrix_approx_eq C_ref C
    @test @test_matrix_approx_eq Ψ_ref Ψ
    @test @test_matrix_approx_eq Π_ref Π
end

### Measurement equation
expect = Dict{Symbol, Matrix}()
h5 = h5open("$path/measurement.h5")
expect[:ZZ] = read(h5, "ZZ")
expect[:DD] = reshape(read(h5, "DD"), 12, 1)
expect[:QQ] = read(h5, "QQ")
expect[:EE] = read(h5, "EE")
close(h5)

TTT, RRR, CCC = solve(model)
actual = measurement(model, TTT, RRR, CCC)
@testset "Checking measurement equation" begin
    for d in (:ZZ, :DD, :QQ, :EE)
        @test @test_matrix_approx_eq expect[d] actual[d]
    end
end

### Pseudo-measurement equation
expect = Dict{Symbol, Any}()
jld = JLD2.jldopen("$path/pseudo_measurement.jld2")
expect[:ZZ_pseudo] = read(jld, "ZZ_pseudo")
expect[:DD_pseudo] = read(jld, "DD_pseudo")
close(jld)

actual = pseudo_measurement(model, TTT, RRR, CCC)
@testset "Checking pseudo-measurement equation" begin
    for d in (:ZZ_pseudo, :DD_pseudo)
        @test @test_matrix_approx_eq expect[d] getfield(actual,d)
    end
end

### Custom settings
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks => Setting(:n_anticipated_shocks, 6))
model = Model990(sub, custom_settings = custom_settings)
@testset "Checking n_ant_shocks == 6" begin
    @test get_setting(model, :n_anticipated_shocks) == 6
end

# Indices initialized correctly under custom settings

@testset "Checking model indices part 2" begin
    # Endogenous states
    endo = model.endogenous_states
    @test length(endo) == 66
    @test endo[:rm_tl6] == 66

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 22
    @test exo[:rm_shl6] == 22

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 13
    @test ex[:Erk_f_sh] == 13

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 66
    @test eq[:eq_rml6] == 66

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 12
    @test endo_new[:y_t1] == 67

    # Observables
    obs = model.observables
    @test length(obs) == 18
    @test obs[:obs_nominalrate6] == 18
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

@testset "Check eqcond and state-space dimensions part 2" begin
    # Matrices are of expected dimensions
    @test size(Γ0) == (66, 66)
    @test size(Γ1) == (66, 66)
    @test size(C) == (66,)
    @test size(Ψ) == (66, 22)
    @test size(Π) == (66, 13)
end

nothing
