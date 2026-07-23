using DSGE, Test, ModelConstructors, BenchmarkTools

# smooth_ait_gdp_alt adds pgap_t and ygap_t and references R_t, π_t, y_t, z_t.
# SmetsWouters is the lightest model with those states and the regime_switching
# augment_states that the solver uses.
m = SmetsWouters()

# Historical-rule conditions and sizes, before the rule augments the states.
Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
n_eq_hist = n_equilibrium_conditions(m)
n_st_hist = n_states(m)
old_eq_mp = m.equilibrium_conditions[:eq_mp]
eq_other  = setdiff(1:n_eq_hist, [old_eq_mp])

Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.smooth_ait_gdp_alt_eqcond(m)

endo = m.endogenous_states
eq   = m.equilibrium_conditions

# Default coefficients (ss0, no overrides)
ρ_pgap   = exp(log(0.5) / 10.)
ρ_ygap   = exp(log(0.5) / 10.)
ρ_smooth = 0.656
φ_π      = 11.13
φ_y      = 11.13

@testset "smooth_ait_gdp_alt AltPolicy object is wired up correctly" begin
    pol = DSGE.smooth_ait_gdp_alt()
    @test pol isa DSGE.AltPolicy
    @test pol.key           == :smooth_ait_gdp_alt
    @test pol.eqcond        === DSGE.smooth_ait_gdp_alt_eqcond
    @test pol.solve         === DSGE.smooth_ait_gdp_alt_solve
    @test pol.forecast_init === DSGE.smooth_ait_gdp_alt_forecast_init
end

@testset "smooth_ait_gdp_alt adds the pgap and ygap states and conditions" begin
    @test haskey(endo, :pgap_t) && haskey(endo, :ygap_t)
    @test haskey(eq, :eq_pgap)  && haskey(eq, :eq_ygap)
    @test endo[:pgap_t] == n_st_hist + 1
    @test endo[:ygap_t] == n_st_hist + 2
    @test eq[:eq_pgap]  == n_eq_hist + 1
    @test eq[:eq_ygap]  == n_eq_hist + 2
    @test size(Γ0_alt, 1) == size(Γ0_alt, 2) == n_states(m) == n_st_hist + 2
    @test size(Γ1_alt) == size(Γ0_alt)
end

@testset "smooth_ait_gdp_alt price-gap law of motion" begin
    # pgap_t = π_t + ρ_pgap·pgap_{t-1}
    @test Γ0_alt[eq[:eq_pgap], endo[:pgap_t]] == 1.
    @test Γ0_alt[eq[:eq_pgap], endo[:π_t]]    == -1.
    @test Γ1_alt[eq[:eq_pgap], endo[:pgap_t]] == ρ_pgap
end

@testset "smooth_ait_gdp_alt GDP-gap law of motion" begin
    @test Γ0_alt[eq[:eq_ygap], endo[:ygap_t]] == 1.
    @test Γ0_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
    @test Γ0_alt[eq[:eq_ygap], endo[:z_t]]    == -1.
    @test Γ1_alt[eq[:eq_ygap], endo[:ygap_t]] == ρ_ygap
    @test Γ1_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
end

@testset "smooth_ait_gdp_alt smoothed MP rule" begin
    @test Γ0_alt[eq[:eq_mp], endo[:R_t]]    == 1.
    @test Γ1_alt[eq[:eq_mp], endo[:R_t]]    == ρ_smooth
    @test Γ0_alt[eq[:eq_mp], endo[:pgap_t]] == -φ_π * (1. - ρ_pgap) * (1. - ρ_smooth)
    @test Γ0_alt[eq[:eq_mp], endo[:ygap_t]] == -φ_y * (1. - ρ_ygap) * (1. - ρ_smooth)
    @test C_alt[eq[:eq_mp]]                 == 0.
    # MP shock is fully removed (no rm_t term)
    @test all(Ψ_alt[eq[:eq_mp], :] .== 0.)
end

@testset "smooth_ait_gdp_alt leaves the non-policy equilibrium conditions unchanged" begin
    @test Γ0_hist[eq_other, 1:n_st_hist] ≈ Γ0_alt[eq_other, 1:n_st_hist]
    @test Γ1_hist[eq_other, 1:n_st_hist] ≈ Γ1_alt[eq_other, 1:n_st_hist]
    @test Ψ_hist[eq_other, :]            ≈ Ψ_alt[eq_other, :]
    @test Π_hist[eq_other, :]            ≈ Π_alt[eq_other, :]
    @test !(Γ0_hist[old_eq_mp, 1:n_st_hist] ≈ Γ0_alt[old_eq_mp, 1:n_st_hist])
end

@testset "ait_Thalf / gdp_Thalf settings feed through to ρ_pgap / ρ_ygap" begin
    mt = SmetsWouters(custom_settings = [Setting(:ait_Thalf, 20.),
                                         Setting(:gdp_Thalf, 5.)])
    Γ0_t, Γ1_t, ~ = DSGE.smooth_ait_gdp_alt_eqcond(mt)
    eqt, endot = mt.equilibrium_conditions, mt.endogenous_states
    @test Γ1_t[eqt[:eq_pgap], endot[:pgap_t]] == exp(log(0.5) / 20.)
    @test Γ1_t[eqt[:eq_ygap], endot[:ygap_t]] == exp(log(0.5) / 5.)
end

@testset "smooth_ait_gdp_alt_solve returns consistently sized, finite transition matrices" begin
    TTT, RRR, CCC = DSGE.smooth_ait_gdp_alt_solve(m)
    n = n_states_augmented(m)
    @test size(TTT) == (n, n)
    @test size(RRR) == (n, n_shocks_exogenous(m))
    @test length(CCC) == n
    @test all(isfinite, TTT)
    @test all(isfinite, RRR)
end

@testset "smooth_ait_gdp_alt_forecast_init seeds pgap / ygap" begin
    m <= Setting(:pgap_value, 2.5)
    m <= Setting(:ygap_value, 1.5)

    n           = n_states_augmented(m)
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n)

    new_shocks, fs = DSGE.smooth_ait_gdp_alt_forecast_init(m, shocks, copy(final_state))
    @test new_shocks == shocks
    @test fs[endo[:pgap_t]] == -2.5
    @test fs[endo[:ygap_t]] == -1.5
    other = setdiff(1:n, [endo[:pgap_t], endo[:ygap_t]])
    @test fs[other] == final_state[other]
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    mb = SmetsWouters()
    mb <= Setting(:pgap_value, 2.5)
    mb <= Setting(:ygap_value, 1.5)

    b_construct = @benchmark DSGE.smooth_ait_gdp_alt()
    b_eqcond    = @benchmark DSGE.smooth_ait_gdp_alt_eqcond($mb)
    b_solve     = @benchmark DSGE.smooth_ait_gdp_alt_solve($mb)

    shocks      = zeros(n_shocks_exogenous(mb), 1)
    final_state = collect(1.0:n_states_augmented(mb))
    b_finit     = @benchmark DSGE.smooth_ait_gdp_alt_forecast_init($mb, $shocks, copy($final_state))

    println("\n===== smooth_ait_gdp_alt benchmark results =====")
    for (name, b) in [("smooth_ait_gdp_alt",               b_construct),
                      ("smooth_ait_gdp_alt_eqcond",        b_eqcond),
                      ("smooth_ait_gdp_alt_solve",         b_solve),
                      ("smooth_ait_gdp_alt_forecast_init", b_finit)]
        println(rpad(name, 36), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
