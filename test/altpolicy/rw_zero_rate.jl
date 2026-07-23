using DSGE, Test, ModelConstructors, BenchmarkTools

m = SmetsWouters()

Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
n_eq_hist = n_equilibrium_conditions(m)
n_st_hist = n_states(m)
old_eq_mp = m.equilibrium_conditions[:eq_mp]
eq_other  = setdiff(1:n_eq_hist, [old_eq_mp])   # non-policy rows

# Build the rw_zero_rate equilibrium conditions.
Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.rw_zero_rate_eqcond(m)

endo = m.endogenous_states
eq   = m.equilibrium_conditions

# Default coefficients
ρ_pgap   = exp(log(0.5) / 10.)  
ρ_ygap   = exp(log(0.5) / 10.)  
ρ_smooth = 0.656
φ_π      = 11.13
φ_y      = 11.13
ρ_rw     = 0.93

@testset "rw_zero_rate AltPolicy object is wired up correctly" begin
    pol = DSGE.rw_zero_rate()
    @test pol isa DSGE.AltPolicy
    @test pol.key           == :rw_zero_rate
    @test pol.eqcond        === DSGE.rw_zero_rate_eqcond
    @test pol.solve         === DSGE.rw_zero_rate_solve
    @test pol.forecast_init === DSGE.rw_zero_rate_forecast_init
end

@testset "rw_zero_rate adds the four states and equilibrium conditions" begin
    for s in (:pgap_t, :ygap_t, :rw_t, :Rref_t)
        @test haskey(endo, s)
    end
    for e in (:eq_pgap, :eq_ygap, :eq_rw, :eq_Rref)
        @test haskey(eq, e)
    end
    # Appended in order: pgap, ygap, rw, Rref
    @test endo[:pgap_t] == n_st_hist + 1
    @test endo[:ygap_t] == n_st_hist + 2
    @test endo[:rw_t]   == n_st_hist + 3
    @test endo[:Rref_t] == n_st_hist + 4
    @test eq[:eq_pgap]  == n_eq_hist + 1
    @test eq[:eq_ygap]  == n_eq_hist + 2
    @test eq[:eq_rw]    == n_eq_hist + 3
    @test eq[:eq_Rref]  == n_eq_hist + 4
    # Augmented eqcond matrices are square at the new (larger) dimension
    @test size(Γ0_alt, 1) == size(Γ0_alt, 2) == n_states(m) == n_st_hist + 4
    @test size(Γ1_alt) == size(Γ0_alt)
end

@testset "rw_zero_rate price-gap law of motion: pgap_t = π_t + ρ_pgap·pgap_{t-1}" begin
    @test Γ0_alt[eq[:eq_pgap], endo[:pgap_t]] == 1.
    @test Γ0_alt[eq[:eq_pgap], endo[:π_t]]    == -1.
    @test Γ1_alt[eq[:eq_pgap], endo[:pgap_t]] == ρ_pgap
end

@testset "rw_zero_rate GDP-gap law of motion" begin
    @test Γ0_alt[eq[:eq_ygap], endo[:ygap_t]] == 1.
    @test Γ0_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
    @test Γ0_alt[eq[:eq_ygap], endo[:z_t]]    == -1.
    @test Γ1_alt[eq[:eq_ygap], endo[:ygap_t]] == ρ_ygap
    @test Γ1_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
end

@testset "rw_zero_rate reference-rate penalty (rw_t)" begin
    # rw_t = ρ_rw·rw_{t-1} + Rref_{t-1} - R_{t-1}
    @test Γ0_alt[eq[:eq_rw], endo[:rw_t]]   == 1.
    @test Γ1_alt[eq[:eq_rw], endo[:rw_t]]   == ρ_rw
    @test Γ1_alt[eq[:eq_rw], endo[:Rref_t]] == 1.
    @test Γ1_alt[eq[:eq_rw], endo[:R_t]]    == -1.
end

@testset "rw_zero_rate reference-rate evolution (Rref_t)" begin
    @test Γ0_alt[eq[:eq_Rref], endo[:Rref_t]] == 1.
    @test Γ1_alt[eq[:eq_Rref], endo[:Rref_t]] == ρ_smooth
    @test C_alt[eq[:eq_Rref]]                 == 0.
    @test Γ0_alt[eq[:eq_Rref], endo[:pgap_t]] == -φ_π * (1. - ρ_pgap) * (1. - ρ_smooth)
    @test Γ0_alt[eq[:eq_Rref], endo[:ygap_t]] == -φ_y * (1. - ρ_ygap) * (1. - ρ_smooth)
end

@testset "rw_zero_rate pins the actual rate at the ZLB" begin
    # Old MP rule wiped, R_t pinned to a constant (the zero lower bound, in
    # deviations from the steady-state nominal rate Rstarn)
    @test Γ0_alt[eq[:eq_mp], endo[:R_t]] == 1.
    @test C_alt[eq[:eq_mp]]              == -m[:Rstarn]
    @test all(Γ1_alt[eq[:eq_mp], :] .== 0.)
    @test all(Ψ_alt[eq[:eq_mp], :]  .== 0.)
    @test all(Π_alt[eq[:eq_mp], :]  .== 0.)
end

@testset "rw_zero_rate leaves the non-policy equilibrium conditions unchanged" begin
    # Compare on the original (pre-augmentation) rows and columns
    @test Γ0_hist[eq_other, 1:n_st_hist] ≈ Γ0_alt[eq_other, 1:n_st_hist]
    @test Γ1_hist[eq_other, 1:n_st_hist] ≈ Γ1_alt[eq_other, 1:n_st_hist]
    @test Ψ_hist[eq_other, :]            ≈ Ψ_alt[eq_other, :]
    @test Π_hist[eq_other, :]            ≈ Π_alt[eq_other, :]
    # ...but the monetary policy row does change
    @test !(Γ0_hist[old_eq_mp, 1:n_st_hist] ≈ Γ0_alt[old_eq_mp, 1:n_st_hist])
end

@testset "ait_Thalf / gdp_Thalf / ρ_rw settings feed through" begin
    mt = SmetsWouters(custom_settings = [Setting(:ait_Thalf, 20.),
                                         Setting(:gdp_Thalf, 5.),
                                         Setting(:ρ_rw, 0.5)])
    Γ0_t, Γ1_t, ~ = DSGE.rw_zero_rate_eqcond(mt)
    eqt   = mt.equilibrium_conditions
    endot = mt.endogenous_states
    @test Γ1_t[eqt[:eq_pgap], endot[:pgap_t]] == exp(log(0.5) / 20.)
    @test Γ1_t[eqt[:eq_ygap], endot[:ygap_t]] == exp(log(0.5) / 5.)
    @test Γ1_t[eqt[:eq_rw],   endot[:rw_t]]   == 0.5
end

@testset "rw_zero_rate_solve is degenerate under permanent application" begin
    # rw_zero_rate permanently pins R_t at the ZLB, removing the Taylor rule that
    # normally delivers determinacy. Applied *permanently* and solved with plain
    # gensys there is no unique stable solution, so the solver throws. Confirmed
    # for both SmetsWouters and Model1002 — this is the rule's nature, not a model
    # choice. (It is meant for *temporary* ZLB regimes via gensys2; rw_zero_rate_solve
    # is never called as a permanent rule anywhere in DSGE.)
    @test_throws DSGE.GensysError DSGE.rw_zero_rate_solve(m)
end

@testset "rw_zero_rate_forecast_init seeds pgap / ygap / rw / Rref" begin
    m <= Setting(:pgap_value, 2.5)
    m <= Setting(:ygap_value, 1.5)
    m <= Setting(:rw_value,   0.75)
    m <= Setting(:Rref_value, 0.1)

    n           = n_states_augmented(m)
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n)            # distinct values to detect mutation

    new_shocks, fs = DSGE.rw_zero_rate_forecast_init(m, shocks, copy(final_state))
    @test new_shocks == shocks                  # shocks pass through untouched
    @test length(fs) == n                       # length preserved
    @test fs[endo[:pgap_t]] == -2.5
    @test fs[endo[:ygap_t]] == -1.5
    @test fs[endo[:rw_t]]   == -0.75
    @test fs[endo[:Rref_t]] == 0.1              # taken from :Rref_value setting
    # Every other state is left untouched
    seeded = [endo[:pgap_t], endo[:ygap_t], endo[:rw_t], endo[:Rref_t]]
    other  = setdiff(1:n, seeded)
    @test fs[other] == final_state[other]
end

@testset "rw_zero_rate_forecast_init falls back to R_t when :Rref_value unset" begin
    mf = SmetsWouters()
    DSGE.rw_zero_rate_eqcond(mf)                # augment so the states exist
    mf <= Setting(:pgap_value, 0.)
    mf <= Setting(:ygap_value, 0.)
    mf <= Setting(:rw_value,   0.)
    # deliberately no :Rref_value

    n           = n_states_augmented(mf)
    shocks      = zeros(n_shocks_exogenous(mf), 1)
    final_state = collect(1.0:n)
    _, fs = DSGE.rw_zero_rate_forecast_init(mf, shocks, copy(final_state))
    @test fs[mf.endogenous_states[:Rref_t]] == final_state[mf.endogenous_states[:R_t]]
end

################
# Benchmarking #
################
# Set this flag to true to run the rw_zero_rate benchmarks. Off by default so
# the test suite stays fast.
run_benchmarks = false

if run_benchmarks
    mb = SmetsWouters()
    mb <= Setting(:pgap_value, 2.5)
    mb <= Setting(:ygap_value, 1.5)
    mb <= Setting(:rw_value,   0.75)
    mb <= Setting(:Rref_value, 0.1)

    # Constructing the AltPolicy object
    b_construct = @benchmark DSGE.rw_zero_rate()

    # Building the rw_zero_rate equilibrium conditions
    b_eqcond = @benchmark DSGE.rw_zero_rate_eqcond($mb)

    # Note: rw_zero_rate_solve is not benchmarked — permanent application has no
    # unique stable solution (gensys throws on every model), so there is nothing
    # to time.

    # Initializing a forecast under the rw_zero_rate rule
    shocks      = zeros(n_shocks_exogenous(mb), 1)
    final_state = collect(1.0:n_states_augmented(mb))
    b_finit = @benchmark DSGE.rw_zero_rate_forecast_init($mb, $shocks, copy($final_state))

    println("\n===== rw_zero_rate benchmark results =====")
    for (name, b) in [("rw_zero_rate",               b_construct),
                      ("rw_zero_rate_eqcond",        b_eqcond),
                      ("rw_zero_rate_forecast_init", b_finit)]
        println(rpad(name, 28), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
