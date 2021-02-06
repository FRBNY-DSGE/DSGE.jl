using DSGE, ModelConstructors, HDF5, Test, Dates
path = dirname(@__FILE__)
#=
### Model
model = Model1002("ss10")

### Parameters

@testset "Test parameter bounds checking" begin
    for θ in model.parameters
        @test isa(θ, AbstractParameter)
        if !θ.fixed
            (left, right) = θ.valuebounds
            @test left < θ.value < right
        end
    end
end


### Model indices

@testset "Checking model indices" begin
    # Endogenous states
    endo = model.endogenous_states
    @test length(endo) == 68
    @test endo[:Ez_t] == 60

    # Exogenous shocks
    exo = model.exogenous_shocks
    @test length(exo) == 24
    @test exo[:corepce_sh] == 16

    # Expectation shocks
    ex = model.expected_shocks
    @test length(ex) == 13
    @test ex[:EL_f_sh] == 12

    # Equations
    eq = model.equilibrium_conditions
    @test length(eq) == 68
    @test eq[:eq_Ez] == 59

    # Additional states
    endo_new = model.endogenous_states_augmented
    @test length(endo_new) == 16
    @test endo_new[:y_t1] == 69

    # Observables
    obs = model.observables
    @test length(obs) == 19
    @test obs[:obs_tfp] == 12
end

### Equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(model)

# Transition and measurement equations
TTT, RRR, CCC = solve(model)
meas = measurement(model, TTT, RRR, CCC)

### Pseudo-measurement equation
pseudo_meas = pseudo_measurement(model, TTT, RRR, CCC)

@testset "Check eqcond and state-space dimensions" begin
    # Matrices are of expected dimensions
    @test size(Γ0) == (68, 68)
    @test size(Γ1) == (68, 68)
    @test size(C) == (68,)
    @test size(Ψ) == (68, 24)
    @test size(Π) == (68, 13)

    @test size(meas[:ZZ]) == (19,84)
    @test size(meas[:DD]) == (19,)
    @test size(meas[:QQ]) == (24,24)
    @test size(meas[:EE]) == (19,19)

    @test size(TTT) == (84,84)
    @test size(RRR) == (84,24)
    @test size(CCC) == (84,)
end

G0_exp, G1_exp, C_exp, PSI_exp, PIE_exp = h5open("$path/eqcond.h5") do file
    read(file, "G0"), read(file, "G1"), read(file, "C"),
    read(file, "PSI"), read(file, "PIE")
end
TTT_exp, RRR_exp, CCC_exp = h5open("$path/transition.h5") do file
    read(file, "TTT"), read(file, "RRR"), read(file, "CCC")
end
Z_exp, D_exp, Q_exp, E_exp = h5open("$path/measurement.h5") do file
    read(file, "ZZ"), read(file, "DD"), read(file, "QQ"), read(file, "EE")
end
Z_pseudo_exp, D_pseudo_exp = h5open("$path/pseudo_measurement.h5") do file
    read(file, "ZZ_pseudo"), read(file, "DD_pseudo")
end

@testset "Compare eqcond, transition, measurement, and pseudo-measurement matrices against reference matrices" begin
    @test @test_matrix_approx_eq Γ0 G0_exp
    @test @test_matrix_approx_eq Γ1 G1_exp
    @test @test_matrix_approx_eq C C_exp
    @test @test_matrix_approx_eq Ψ PSI_exp
    @test @test_matrix_approx_eq Π PIE_exp

    @test @test_matrix_approx_eq TTT TTT_exp
    @test @test_matrix_approx_eq RRR RRR_exp
    @test @test_matrix_approx_eq CCC CCC_exp

    @test @test_matrix_approx_eq meas[:ZZ] Z_exp
    @test @test_matrix_approx_eq meas[:DD] D_exp
    @test @test_matrix_approx_eq meas[:QQ] Q_exp
    @test @test_matrix_approx_eq meas[:EE] E_exp

    @test @test_matrix_approx_eq pseudo_meas[:ZZ_pseudo] Z_pseudo_exp[1:21, :] # Test files had been generated w/extra pseudo-obs
    @test @test_matrix_approx_eq pseudo_meas[:DD_pseudo] D_pseudo_exp[1:21, :] # Test files had been generated w/extra pseudo-obs
end
=#
# Check ss62 works
m = Model1002("ss62")
for (i, d) in enumerate([date_presample_start(m), Date(2020, 3, 31), Date(2020, 6, 30), Date(2020, 9, 30),
                         Date(2020, 12, 31), Date(2021, 3, 31)])
    @test get_setting(m, :regime_dates)[i] == d
end
sys = compute_system(m)
exo = m.exogenous_shocks
inds1 = 1:24 # standard shock indices
inds2 = vcat(1:14, 17:24) # standard shock indices excluding measurement error
inds3 = 25:30 # covid shock indices
for i in 2:get_setting(m, :n_regimes)
    @test sys[1, :TTT] ≈ sys[i, :TTT]
    # @test sys[1, :RRR] ≈ sys[i, :RRR]
    @test sys[1, :CCC] ≈ sys[i, :CCC]
    @test sys[1, :ZZ]  ≈ sys[i, :ZZ]
    @test sys[1, :DD]  ≈ sys[i, :DD]

    if i >= 5
        @test sys[1, :QQ][inds1, inds1] ≈ sys[i, :QQ][inds1, inds1]
    end
    if i >= 4
        @test sys[1, :QQ][inds2, inds2] ≈ sys[i, :QQ][inds2, inds2]
    end
    if 2 <= i <= 5
        @test sys[i, :QQ][exo[:biidc_sh], exo[:biidc_sh]] == 16
        @test sys[i, :QQ][exo[:ziid_sh], exo[:ziid_sh]] == 25.
        @test sys[i, :QQ][exo[:φ_sh], exo[:φ_sh]] == 160000.
    end
    @test sys[i, :QQ][exo[:ziid_shl1], exo[:ziid_shl1]] == 0.
    @test sys[i, :QQ][exo[:φ_shl1], exo[:φ_shl1]] == 0.
    if 2 <= i <= 3 || i == 6
        @test sys[i, :QQ][exo[:biidc_shl1], exo[:biidc_shl1]] == 0.
    else
        @test sys[i, :QQ][exo[:biidc_shl1], exo[:biidc_shl1]] == 16.
    end
end
