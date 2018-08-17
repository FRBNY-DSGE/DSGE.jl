using DSGE
using HDF5, Base.Test, Distributions

path = dirname(@__FILE__)

### Model
sw = SmetsWouters()

### Parameters

@testset "Test parameter bounds checking" begin
    for θ in sw.parameters
        !isa(θ,AbstractParameter) && error()
        θ.fixed && continue

        (left, right) = θ.valuebounds
        @test left < θ.value < right

    end
end

### Model indices

@testset "Checking model indices" begin
    # Endogenous states
    endo = sw.endogenous_states
    @test length(endo) == 47
    @test endo[:ztil_t] == 47

    # Exogenous shocks
    exo = sw.exogenous_shocks
    @test length(exo) == 7
    @test exo[:rm_sh] == 7

    # Expectation shocks
    ex = sw.expected_shocks
    @test length(ex) == 12
    @test ex[:Erk_f_sh] == 12

    # Equations
    eq = sw.equilibrium_conditions
    @test length(eq) == 47
    @test eq[:eq_ztil] == 47

    # Additional states
    endo_new = sw.endogenous_states_augmented
    @test length(endo_new) == 7
    @test endo_new[:y_t1] == 48

    # Observables
    obs = sw.observables
    @test length(obs) == 7
    @test obs[:obs_investment] == 7
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(sw)

@testset "Check eqcond and state-space dimensions" begin
    # Matrices are of expected dimensions
    @test size(Γ0) == (47, 47)
    @test size(Γ1) == (47, 47)
    @test size(C) == (47,)
    @test size(Ψ) == (47, 7)
    @test size(Π) == (47, 12)
end

# Check output matrices against reference output (ϵ = 1e-4)
h5 = h5open("$path/eqcond.h5")
Γ0_ref = read(h5, "G0")
Γ1_ref = read(h5, "G1")
C_ref  = reshape(read(h5, "C"), 47, 1)
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
expect[:DD] = reshape(read(h5, "DD"), 7, 1)
expect[:QQ] = read(h5, "QQ")
expect[:EE] = read(h5, "EE")
close(h5)

sw = SmetsWouters()
TTT, RRR, CCC = solve(sw)
actual = measurement(sw, TTT, RRR, CCC)
@testset "Checking measurement equation" begin
    for d in (:ZZ, :DD, :QQ, :EE)
        @test @test_matrix_approx_eq expect[d] actual[d]
    end
end

### Custom settings
custom_settings = Dict{Symbol, Setting}(
    :reoptimize => Setting(:reoptimize, false))
model = SmetsWouters(custom_settings = custom_settings)
@testset "Check SW custom setting" begin
    @test get_setting(model, :reoptimize) == false
end

nothing
