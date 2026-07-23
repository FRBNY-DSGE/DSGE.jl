using DSGE, Test, ModelConstructors, BenchmarkTools

# taylor93 replaces the MP rule with Taylor's (1993) rule, which references the
# annualized inflation π_a_t and flexible-price output y_f_t. Those states only
# exist in the full models, so we use Model990. No states are added.
m = Model990()

Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
eq, endo = m.equilibrium_conditions, m.endogenous_states
eq_mp    = eq[:eq_mp]
eq_other = setdiff(1:n_equilibrium_conditions(m), [eq_mp])

Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.taylor93_eqcond(m)

@testset "taylor93 AltPolicy object is wired up correctly" begin
    pol = DSGE.taylor93()
    @test pol.key           == :taylor93
    @test pol.eqcond        === DSGE.taylor93_eqcond
    @test pol.solve         === DSGE.taylor93_solve
    @test pol.forecast_init === identity   # no forecast_init given
end

@testset "taylor93 keeps the state space the same size" begin
    @test size(Γ0_alt) == size(Γ0_hist)
    @test size(Γ1_alt) == size(Γ1_hist)
end

@testset "taylor93 rule has the right coefficients" begin
    @test Γ0_alt[eq_mp, endo[:R_t]]   == 1
    @test Γ0_alt[eq_mp, endo[:π_a_t]] == -1.5/4
    @test Γ0_alt[eq_mp, endo[:y_t]]   == -0.5/4
    @test Γ0_alt[eq_mp, endo[:y_f_t]] == 0.5/4
    # Old rule fully zeroed out
    @test all(Γ1_alt[eq_mp, :] .== 0)
    @test C_alt[eq_mp] == 0
    @test all(Ψ_alt[eq_mp, :] .== 0)
    @test all(Π_alt[eq_mp, :] .== 0)
end

@testset "taylor93 leaves the non-policy equilibrium conditions unchanged" begin
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt[eq_other, :]
    @test !(Γ0_hist[eq_mp, :] ≈ Γ0_alt[eq_mp, :])
end

@testset "taylor93_solve returns consistently sized, finite transition matrices" begin
    TTT, RRR, CCC = DSGE.taylor93_solve(m)
    n = n_states_augmented(m)
    @test size(TTT) == (n, n)
    @test size(RRR) == (n, n_shocks_exogenous(m))
    @test length(CCC) == n
    @test all(isfinite, TTT)
    @test all(isfinite, RRR)
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    mb = Model990()
    b_construct = @benchmark DSGE.taylor93()
    b_eqcond    = @benchmark DSGE.taylor93_eqcond($mb)
    b_solve     = @benchmark DSGE.taylor93_solve($mb)

    println("\n===== taylor93 benchmark results =====")
    for (name, b) in [("taylor93",        b_construct),
                      ("taylor93_eqcond", b_eqcond),
                      ("taylor93_solve",  b_solve)]
        println(rpad(name, 18), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
