using DSGE, Test, ModelConstructors, BenchmarkTools

m = SmetsWouters()
Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
n_eq_hist = n_equilibrium_conditions(m)
n_st_hist = n_states(m)
old_eq_mp = m.equilibrium_conditions[:eq_mp]
eq_other  = setdiff(1:n_eq_hist, [old_eq_mp])   # non-policy rows


Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.flexible_ait_eqcond(m)

endo = m.endogenous_states
eq   = m.equilibrium_conditions

ρ_pgap   = exp(log(0.5) / 10.)   # ait_Thalf defaults to 10
ρ_ygap   = exp(log(0.5) / 10.)   # gdp_Thalf defaults to 10
ρ_smooth = 0.656
φ_π      = 11.13
φ_y      = 11.13

@testset "flexible_ait AltPolicy object is wired up correctly" begin
    pol = DSGE.flexible_ait()
    @test pol isa DSGE.AltPolicy
    @test pol.key           == :flexible_ait
    @test pol.eqcond        === DSGE.flexible_ait_eqcond
    @test pol.solve         === DSGE.flexible_ait_solve
    @test pol.forecast_init === DSGE.flexible_ait_forecast_init
end

@testset "flexible_ait adds the pgap and ygap states and conditions" begin
    @test haskey(endo, :pgap_t) && haskey(endo, :ygap_t)
    @test haskey(eq, :eq_pgap)  && haskey(eq, :eq_ygap)
    @test endo[:pgap_t] == n_st_hist + 1
    @test endo[:ygap_t] == n_st_hist + 2
    @test eq[:eq_pgap]  == n_eq_hist + 1
    @test eq[:eq_ygap]  == n_eq_hist + 2
    # Augmented eqcond matrices are square at the new (larger) dimension
    @test size(Γ0_alt, 1) == size(Γ0_alt, 2) == n_states(m) == n_st_hist + 2
    @test size(Γ1_alt) == size(Γ0_alt)
end

@testset "flexible_ait price-gap law of motion: pgap_t = π_t + ρ_pgap·pgap_{t-1}" begin
    @test Γ0_alt[eq[:eq_pgap], endo[:pgap_t]] == 1.
    @test Γ0_alt[eq[:eq_pgap], endo[:π_t]]    == -1.
    @test Γ1_alt[eq[:eq_pgap], endo[:pgap_t]] == ρ_pgap
end

@testset "flexible_ait GDP-gap law of motion" begin
    @test Γ0_alt[eq[:eq_ygap], endo[:ygap_t]] == 1.
    @test Γ0_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
    @test Γ0_alt[eq[:eq_ygap], endo[:z_t]]    == -1.
    @test Γ1_alt[eq[:eq_ygap], endo[:ygap_t]] == ρ_ygap
    @test Γ1_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
end

@testset "flexible_ait monetary policy rule has the right coefficients" begin
    @test Γ0_alt[eq[:eq_mp], endo[:R_t]]    == 1.
    @test Γ1_alt[eq[:eq_mp], endo[:R_t]]    == ρ_smooth
    @test Γ0_alt[eq[:eq_mp], endo[:pgap_t]] == -φ_π * (1. - ρ_pgap) * (1. - ρ_smooth)
    @test Γ0_alt[eq[:eq_mp], endo[:ygap_t]] == -φ_y * (1. - ρ_ygap) * (1. - ρ_smooth)
    @test C_alt[eq[:eq_mp]]                 == 0.
    # No add_ait_rm → the monetary policy shock enters via rm_t
    @test Γ0_alt[eq[:eq_mp], endo[:rm_t]]   == -1.
end

@testset "flexible_ait leaves the non-policy equilibrium conditions unchanged" begin
    # Compare on the original (pre-augmentation) rows and columns
    @test Γ0_hist[eq_other, 1:n_st_hist] ≈ Γ0_alt[eq_other, 1:n_st_hist]
    @test Γ1_hist[eq_other, 1:n_st_hist] ≈ Γ1_alt[eq_other, 1:n_st_hist]
    @test Ψ_hist[eq_other, :]            ≈ Ψ_alt[eq_other, :]
    @test Π_hist[eq_other, :]            ≈ Π_alt[eq_other, :]
    # ...but the monetary policy row does change
    @test !(Γ0_hist[old_eq_mp, 1:n_st_hist] ≈ Γ0_alt[old_eq_mp, 1:n_st_hist])
end

@testset "ait_Thalf / gdp_Thalf settings feed through to ρ_pgap / ρ_ygap" begin
    mt = SmetsWouters(custom_settings = [Setting(:ait_Thalf, 20.),
                                         Setting(:gdp_Thalf, 5.)])
    Γ0_t, Γ1_t, ~ = DSGE.flexible_ait_eqcond(mt)
    @test Γ1_t[mt.equilibrium_conditions[:eq_pgap], mt.endogenous_states[:pgap_t]] ==
        exp(log(0.5) / 20.)
    @test Γ1_t[mt.equilibrium_conditions[:eq_ygap], mt.endogenous_states[:ygap_t]] ==
        exp(log(0.5) / 5.)
end

@testset "flexible_ait_solve returns consistently sized, finite transition matrices" begin
    TTT, RRR, CCC = DSGE.flexible_ait_solve(m)
    n = n_states_augmented(m)
    @test size(TTT) == (n, n)
    @test size(RRR) == (n, n_shocks_exogenous(m))
    @test length(CCC) == n
    @test all(isfinite, TTT)
    @test all(isfinite, RRR)
end

@testset "flexible_ait_forecast_init sets pgap / ygap to their target values" begin
    m <= Setting(:pgap_value, 2.5)
    m <= Setting(:ygap_value, 1.5)

    n           = n_states_augmented(m)
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n)            # distinct values to detect mutation

    new_shocks, fs = DSGE.flexible_ait_forecast_init(m, shocks, copy(final_state))
    @test new_shocks == shocks               # shocks pass through untouched
    @test fs[endo[:pgap_t]] == -2.5
    @test fs[endo[:ygap_t]] == -1.5
    # Every other state is left untouched
    other = setdiff(1:n, [endo[:pgap_t], endo[:ygap_t]])
    @test fs[other] == final_state[other]
end

################
# Benchmarking #
################
# Set this flag to true to run the flexible_ait benchmarks. Off by default so
# the test suite stays fast.
run_benchmarks = false

if run_benchmarks
    mb = SmetsWouters()
    mb <= Setting(:pgap_value, 0.)
    mb <= Setting(:ygap_value, 0.)

    # Constructing the AltPolicy object
    b_construct = @benchmark DSGE.flexible_ait()

    # Building the flexible_ait equilibrium conditions
    b_eqcond = @benchmark DSGE.flexible_ait_eqcond($mb)

    # Solving the model under the flexible_ait rule (eqcond + gensys + augment_states)
    b_solve = @benchmark DSGE.flexible_ait_solve($mb)

    # Initializing a forecast under the flexible_ait rule
    shocks      = zeros(n_shocks_exogenous(mb), 1)
    final_state = collect(1.0:n_states_augmented(mb))
    b_finit = @benchmark DSGE.flexible_ait_forecast_init($mb, $shocks, copy($final_state))

    println("\n===== flexible_ait benchmark results =====")
    for (name, b) in [("flexible_ait",               b_construct),
                      ("flexible_ait_eqcond",        b_eqcond),
                      ("flexible_ait_solve",         b_solve),
                      ("flexible_ait_forecast_init", b_finit)]
        println(rpad(name, 28), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
