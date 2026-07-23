using DSGE, Test, ModelConstructors, BenchmarkTools

# taylor_rule replaces the MP rule with a smoothed Taylor rule whose coefficients
# come from the parameters ρ, ψ1, ψ2, ψ3. It references π_star_t, y_f_t and the
# MP-shock state rm_t, which only the full models have, so we use Model990. No
# states are added.
m = Model990()

Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
eq, endo = m.equilibrium_conditions, m.endogenous_states
eq_mp    = eq[:eq_mp]
eq_other = setdiff(1:n_equilibrium_conditions(m), [eq_mp])

Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.taylor_rule_eqcond(m)

ρ, ψ1, ψ2, ψ3 = m[:ρ].value, m[:ψ1].value, m[:ψ2].value, m[:ψ3].value

@testset "taylor_rule AltPolicy object is wired up correctly" begin
    pol = DSGE.taylor_rule()
    @test pol.key           == :taylor_rule
    @test pol.eqcond        === DSGE.taylor_rule_eqcond
    @test pol.solve         === DSGE.taylor_rule_solve
    @test pol.forecast_init === DSGE.taylor_rule_forecast_init
end

@testset "taylor_rule keeps the state space the same size" begin
    @test size(Γ0_alt) == size(Γ0_hist)
    @test size(Γ1_alt) == size(Γ1_hist)
end

@testset "taylor_rule has the right coefficients" begin
    @test Γ0_alt[eq_mp, endo[:R_t]]      == 1.
    @test Γ1_alt[eq_mp, endo[:R_t]]      == ρ
    @test C_alt[eq_mp]                   == 0.
    @test Γ0_alt[eq_mp, endo[:π_t]]      == -(1. - ρ) * ψ1
    @test Γ0_alt[eq_mp, endo[:π_star_t]] == (1. - ρ) * ψ1
    @test Γ0_alt[eq_mp, endo[:y_t]]      == -(1. - ρ) * ψ2 - ψ3
    @test Γ0_alt[eq_mp, endo[:y_f_t]]    == (1. - ρ) * ψ2 + ψ3
    @test Γ1_alt[eq_mp, endo[:y_t]]      == -ψ3
    @test Γ1_alt[eq_mp, endo[:y_f_t]]    == ψ3
    @test Γ0_alt[eq_mp, endo[:rm_t]]     == -1.   # MP shock retained
end

@testset "taylor_rule reproduces the model's baseline Taylor rule" begin
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt[eq_other, :]
    @test Γ0_hist[eq_mp, :] ≈ Γ0_alt[eq_mp, :]
    @test Γ1_hist[eq_mp, :] ≈ Γ1_alt[eq_mp, :]
end

@testset "taylor_rule_solve returns consistently sized, finite transition matrices" begin
    TTT, RRR, CCC = DSGE.taylor_rule_solve(m)
    n = n_states_augmented(m)
    @test size(TTT) == (n, n)
    @test size(RRR) == (n, n_shocks_exogenous(m))
    @test length(CCC) == n
    @test all(isfinite, TTT)
    @test all(isfinite, RRR)
end

@testset "taylor_rule_forecast_init passes shocks and state through unchanged" begin
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n_states_augmented(m))
    new_shocks, fs = DSGE.taylor_rule_forecast_init(m, shocks, copy(final_state))
    @test new_shocks == shocks
    @test fs == final_state
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    mb = Model990()
    b_construct = @benchmark DSGE.taylor_rule()
    b_eqcond    = @benchmark DSGE.taylor_rule_eqcond($mb)
    b_solve     = @benchmark DSGE.taylor_rule_solve($mb)

    shocks      = zeros(n_shocks_exogenous(mb), 1)
    final_state = collect(1.0:n_states_augmented(mb))
    b_finit     = @benchmark DSGE.taylor_rule_forecast_init($mb, $shocks, copy($final_state))

    println("\n===== taylor_rule benchmark results =====")
    for (name, b) in [("taylor_rule",               b_construct),
                      ("taylor_rule_eqcond",        b_eqcond),
                      ("taylor_rule_solve",         b_solve),
                      ("taylor_rule_forecast_init", b_finit)]
        println(rpad(name, 28), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
