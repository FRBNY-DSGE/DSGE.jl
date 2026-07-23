using DSGE, Test, ModelConstructors, BenchmarkTools


m = Model990()

Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)

eq   = m.equilibrium_conditions
endo = m.endogenous_states
eq_mp    = eq[:eq_mp]
eq_other = setdiff(1:n_equilibrium_conditions(m), [eq_mp])

Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.alt_inflation_eqcond(m)

@testset "alt_inflation keeps the state space the same size" begin
    @test size(Γ0_alt) == size(Γ0_hist)
    @test size(Γ1_alt) == size(Γ1_hist)
    @test length(C_alt) == length(C_hist)
end

@testset "alt_inflation monetary policy rule: π_t = π_star_t" begin
    # New rule: π_t - π_star_t = 0
    @test Γ0_alt[eq_mp, endo[:π_t]]      == 1
    @test Γ0_alt[eq_mp, endo[:π_star_t]] == -1
    # Every other entry of the MP row is zeroed out
    others = setdiff(1:n_states(m), [endo[:π_t], endo[:π_star_t]])
    @test all(Γ0_alt[eq_mp, others] .== 0)
    @test all(Γ1_alt[eq_mp, :] .== 0)
    @test C_alt[eq_mp] == 0
    @test all(Ψ_alt[eq_mp, :] .== 0)
    @test all(Π_alt[eq_mp, :] .== 0)
end

@testset "alt_inflation leaves the non-policy equilibrium conditions unchanged" begin
    @test Γ0_hist[eq_other, :] ≈ Γ0_alt[eq_other, :]
    @test Γ1_hist[eq_other, :] ≈ Γ1_alt[eq_other, :]
    @test Ψ_hist[eq_other, :] ≈ Ψ_alt[eq_other, :]
    @test Π_hist[eq_other, :] ≈ Π_alt[eq_other, :]
    # ...but the monetary policy row does change
    @test !(Γ0_hist[eq_mp, :] ≈ Γ0_alt[eq_mp, :])
end

@testset "alt_inflation_solve returns consistently sized transition matrices" begin
    TTT, RRR, CCC = DSGE.alt_inflation_solve(m)
    n = n_states_augmented(m)
    @test size(TTT) == (n, n)
    @test size(RRR) == (n, n_shocks_exogenous(m))
    @test length(CCC) == n
    @test all(isfinite, TTT)
    @test all(isfinite, RRR)
end

@testset "alt_inflation AltPolicy object is wired up correctly" begin
    pol = DSGE.alt_inflation()
    @test pol.key == :alt_inflation
    @test pol.eqcond === DSGE.alt_inflation_eqcond
    @test pol.solve === DSGE.alt_inflation_solve
end

################
# Benchmarking #
################
# Set this flag to true to run the alt_inflation benchmarks. Off by default so
# the test suite stays fast.
run_benchmarks = false

if run_benchmarks
    mb = Model990()

    # Building the alt_inflation equilibrium conditions
    b_eqcond = @benchmark DSGE.alt_inflation_eqcond($mb)

    # Solving the model under the alt_inflation rule (eqcond + gensys + augment_states)
    b_solve = @benchmark DSGE.alt_inflation_solve($mb)

    println("\n===== alt_inflation benchmark results =====")
    for (name, b) in [("alt_inflation_eqcond", b_eqcond),
                      ("alt_inflation_solve", b_solve)]
        println(rpad(name, 22), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
