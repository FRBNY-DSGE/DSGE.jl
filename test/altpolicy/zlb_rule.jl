using DSGE, Test, ModelConstructors, BenchmarkTools

# zlb_rule pins R_t at a constant (zlb value / 4 - Rstarn)
m = SmetsWouters()
m <= Setting(:zlb_rule_remove_mon_anticipated_shocks, false)
m <= Setting(:zlb_rule_value, 0.2)        # known value so we can check C

Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
eq, endo = m.equilibrium_conditions, m.endogenous_states
eq_mp    = eq[:eq_mp]
eq_other = setdiff(1:n_equilibrium_conditions(m), [eq_mp])

Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.zlb_rule_eqcond(m)

@testset "zlb_rule AltPolicy object is wired up correctly" begin
    pol = DSGE.zlb_rule()
    @test pol.key           == :zlb_rule
    @test pol.eqcond        === DSGE.zlb_rule_eqcond
    @test pol.solve         === DSGE.zlb_rule_solve
    @test pol.forecast_init === DSGE.zlb_rule_forecast_init
end

@testset "zlb_rule pins R_t at the ZLB constant" begin
    @test size(Γ0_alt) == size(Γ0_hist)
    @test Γ0_alt[eq_mp, endo[:R_t]] == 1.
    @test C_alt[eq_mp]              == 0.2/4 - m[:Rstarn]
    @test all(Γ1_alt[eq_mp, :] .== 0.)
end

@testset "zlb_rule leaves the non-policy equilibrium conditions unchanged" begin
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt[eq_other, :]
    @test !(Γ0_hist[eq_mp, :] ≈ Γ0_alt[eq_mp, :])
end

@testset "zlb_rule_solve is degenerate (interest-rate peg)" begin
    # Pinning R_t to a constant removes the Taylor rule, so there is no unique
    # stable solution and gensys throws. This rule is meant for temporary ZLB
    # regimes via gensys2, not standalone solving.
    @test_throws DSGE.GensysError DSGE.zlb_rule_solve(m)
end

@testset "zlb_rule_forecast_init passes shocks and state through unchanged" begin
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n_states_augmented(m))
    new_shocks, fs = DSGE.zlb_rule_forecast_init(m, shocks, copy(final_state))
    @test new_shocks == shocks
    @test fs == final_state
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    mb = SmetsWouters()
    mb <= Setting(:zlb_rule_remove_mon_anticipated_shocks, false)
    mb <= Setting(:zlb_rule_value, 0.2)
    b_construct = @benchmark DSGE.zlb_rule()
    b_eqcond    = @benchmark DSGE.zlb_rule_eqcond($mb)
    # zlb_rule_solve is not benchmarked — the peg has no stable solution (throws).

    shocks      = zeros(n_shocks_exogenous(mb), 1)
    final_state = collect(1.0:n_states_augmented(mb))
    b_finit     = @benchmark DSGE.zlb_rule_forecast_init($mb, $shocks, copy($final_state))

    println("\n===== zlb_rule benchmark results =====")
    for (name, b) in [("zlb_rule",               b_construct),
                      ("zlb_rule_eqcond",        b_eqcond),
                      ("zlb_rule_forecast_init", b_finit)]
        println(rpad(name, 24), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
