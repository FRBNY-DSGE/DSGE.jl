using DSGE, Test, ModelConstructors, BenchmarkTools

# ngdp (nominal GDP targeting) adds a single price-gap state and references only
# R_t, π_t, y_t and z_t in its rule. We use SmetsWouters rather than the lighter
# AnSchorfheide because ngdp_solve calls augment_states with the regime_switching
# keyword, which only the larger models (SmetsWouters, Model990, Model1002, ...)
# accept. SmetsWouters is the lightest model that supports it.
m = SmetsWouters()

Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
n_eq_hist = n_equilibrium_conditions(m)
n_st_hist = n_states(m)
old_eq_mp = m.equilibrium_conditions[:eq_mp]
eq_other  = setdiff(1:n_eq_hist, [old_eq_mp])   # non-policy rows

Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.ngdp_eqcond(m)

endo = m.endogenous_states
eq   = m.equilibrium_conditions

# Hardcoded ngdp coefficients (see ngdp_replace_eq_entries)
ρ  = 0.
φ1 = 0.25

@testset "ngdp AltPolicy object is wired up correctly" begin
    pol = DSGE.ngdp()
    @test pol isa DSGE.AltPolicy
    @test pol.key           == :ngdp
    @test pol.eqcond        === DSGE.ngdp_eqcond
    @test pol.solve         === DSGE.ngdp_solve
    @test pol.forecast_init === DSGE.ngdp_forecast_init
end

@testset "ngdp adds the pgap state and equilibrium condition" begin
    @test haskey(endo, :pgap_t)
    @test haskey(eq, :eq_pgap)
    @test endo[:pgap_t] == n_st_hist + 1
    @test eq[:eq_pgap]  == n_eq_hist + 1
    # Augmented eqcond matrices are square at the new (larger) dimension
    @test size(Γ0_alt, 1) == size(Γ0_alt, 2) == n_states(m) == n_st_hist + 1
    @test size(Γ1_alt) == size(Γ0_alt)
end

@testset "ngdp price-gap (nominal GDP) law of motion" begin
    # pgap_t = pgap_{t-1} + π_t + y_t + z_t - y_{t-1}
    @test Γ0_alt[eq[:eq_pgap], endo[:pgap_t]] == 1.
    @test Γ0_alt[eq[:eq_pgap], endo[:π_t]]    == -1.
    @test Γ1_alt[eq[:eq_pgap], endo[:pgap_t]] == 1.
    @test Γ0_alt[eq[:eq_pgap], endo[:y_t]]    == -1.
    @test Γ0_alt[eq[:eq_pgap], endo[:z_t]]    == -1.
    @test Γ1_alt[eq[:eq_pgap], endo[:y_t]]    == -1.
end

@testset "ngdp monetary policy rule has the right coefficients" begin
    @test Γ0_alt[eq[:eq_mp], endo[:R_t]]    == 1.
    @test Γ0_alt[eq[:eq_mp], endo[:pgap_t]] == -φ1
    @test Γ1_alt[eq[:eq_mp], endo[:R_t]]    == ρ
    @test C_alt[eq[:eq_mp]]                 == 0.
    # The old monetary policy rule's shock loadings are fully zeroed out
    @test all(Ψ_alt[eq[:eq_mp], :] .== 0.)
end

@testset "ngdp leaves the non-policy equilibrium conditions unchanged" begin
    # Compare on the original (pre-augmentation) rows and columns
    @test Γ0_hist[eq_other, 1:n_st_hist] ≈ Γ0_alt[eq_other, 1:n_st_hist]
    @test Γ1_hist[eq_other, 1:n_st_hist] ≈ Γ1_alt[eq_other, 1:n_st_hist]
    @test Ψ_hist[eq_other, :]            ≈ Ψ_alt[eq_other, :]
    @test Π_hist[eq_other, :]            ≈ Π_alt[eq_other, :]
    # ...but the monetary policy row does change
    @test !(Γ0_hist[old_eq_mp, 1:n_st_hist] ≈ Γ0_alt[old_eq_mp, 1:n_st_hist])
end

@testset "ngdp_solve returns consistently sized, finite transition matrices" begin
    TTT, RRR, CCC = DSGE.ngdp_solve(m)
    n = n_states_augmented(m)
    @test size(TTT) == (n, n)
    @test size(RRR) == (n, n_shocks_exogenous(m))
    @test length(CCC) == n
    @test all(isfinite, TTT)
    @test all(isfinite, RRR)
end

@testset "ngdp_forecast_init sets pgap to -pgap_value" begin
    m <= Setting(:pgap_value, 12.)

    n           = n_states_augmented(m)
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n)            # distinct values to detect mutation

    new_shocks, fs = DSGE.ngdp_forecast_init(m, shocks, copy(final_state))
    @test new_shocks == shocks               # shocks pass through untouched
    @test fs[endo[:pgap_t]] == -12.
    # Every other state is left untouched
    other = setdiff(1:n, [endo[:pgap_t]])
    @test fs[other] == final_state[other]
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    mb = SmetsWouters()
    mb <= Setting(:pgap_value, 12.)

    # Constructing the AltPolicy object
    b_construct = @benchmark DSGE.ngdp()

    # Building the ngdp equilibrium conditions
    b_eqcond = @benchmark DSGE.ngdp_eqcond($mb)

    # Solving the model under the ngdp rule (eqcond + gensys + augment_states)
    b_solve = @benchmark DSGE.ngdp_solve($mb)

    # Initializing a forecast under the ngdp rule
    shocks      = zeros(n_shocks_exogenous(mb), 1)
    final_state = collect(1.0:n_states_augmented(mb))
    b_finit = @benchmark DSGE.ngdp_forecast_init($mb, $shocks, copy($final_state))

    println("\n===== ngdp benchmark results =====")
    for (name, b) in [("ngdp",               b_construct),
                      ("ngdp_eqcond",        b_eqcond),
                      ("ngdp_solve",         b_solve),
                      ("ngdp_forecast_init", b_finit)]
        println(rpad(name, 20), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
