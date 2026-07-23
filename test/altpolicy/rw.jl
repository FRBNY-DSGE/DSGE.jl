using DSGE, Test, ModelConstructors, BenchmarkTools


m = SmetsWouters()

# Historical-rule conditions and sizes, recorded before rw augments the states.
Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
n_eq_hist = n_equilibrium_conditions(m)
n_st_hist = n_states(m)
old_eq_mp = m.equilibrium_conditions[:eq_mp]
eq_other  = setdiff(1:n_eq_hist, [old_eq_mp])

Γ0_alt, Γ1_alt, C_alt, Ψ_alt, Π_alt = DSGE.rw_eqcond(m)

endo = m.endogenous_states
eq   = m.equilibrium_conditions

# Default coefficients (ss0, no overrides)
ρ_pgap   = exp(log(0.5) / 10.)
ρ_ygap   = exp(log(0.5) / 10.)
ρ_smooth = 0.656
φ_π      = 11.13
φ_y      = 11.13
ρ_rw     = 0.93

@testset "rw AltPolicy object is wired up correctly" begin
    pol = DSGE.rw()
    @test pol isa DSGE.AltPolicy
    @test pol.key           == :rw
    @test pol.eqcond        === DSGE.rw_eqcond
    @test pol.solve         === DSGE.rw_solve
    @test pol.forecast_init === DSGE.rw_forecast_init
end

@testset "rw adds the four states and equilibrium conditions" begin
    for s in (:pgap_t, :ygap_t, :rw_t, :Rref_t); @test haskey(endo, s); end
    for e in (:eq_pgap, :eq_ygap, :eq_rw, :eq_Rref); @test haskey(eq, e); end
    # Appended in order: pgap, ygap, rw, Rref
    @test endo[:pgap_t] == n_st_hist + 1
    @test endo[:ygap_t] == n_st_hist + 2
    @test endo[:rw_t]   == n_st_hist + 3
    @test endo[:Rref_t] == n_st_hist + 4
    @test eq[:eq_pgap]  == n_eq_hist + 1
    @test eq[:eq_ygap]  == n_eq_hist + 2
    @test eq[:eq_rw]    == n_eq_hist + 3
    @test eq[:eq_Rref]  == n_eq_hist + 4
    @test size(Γ0_alt, 1) == size(Γ0_alt, 2) == n_states(m) == n_st_hist + 4
    @test size(Γ1_alt) == size(Γ0_alt)
end

@testset "rw price-gap law of motion" begin
    # pgap_t = π_t + ρ_pgap·pgap_{t-1}
    @test Γ0_alt[eq[:eq_pgap], endo[:pgap_t]] == 1.
    @test Γ0_alt[eq[:eq_pgap], endo[:π_t]]    == -1.
    @test Γ1_alt[eq[:eq_pgap], endo[:pgap_t]] == ρ_pgap
end

@testset "rw GDP-gap law of motion" begin
    @test Γ0_alt[eq[:eq_ygap], endo[:ygap_t]] == 1.
    @test Γ0_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
    @test Γ0_alt[eq[:eq_ygap], endo[:z_t]]    == -1.
    @test Γ1_alt[eq[:eq_ygap], endo[:ygap_t]] == ρ_ygap
    @test Γ1_alt[eq[:eq_ygap], endo[:y_t]]    == -1.
end

@testset "rw reference-rate penalty (rw_t)" begin
    # rw_t = ρ_rw·rw_{t-1} + Rref_{t-1} - R_{t-1}
    @test Γ0_alt[eq[:eq_rw], endo[:rw_t]]   == 1.
    @test Γ1_alt[eq[:eq_rw], endo[:rw_t]]   == ρ_rw
    @test Γ1_alt[eq[:eq_rw], endo[:Rref_t]] == 1.
    @test Γ1_alt[eq[:eq_rw], endo[:R_t]]    == -1.
end

@testset "rw reference-rate evolution (Rref_t)" begin
    @test Γ0_alt[eq[:eq_Rref], endo[:Rref_t]] == 1.
    @test Γ1_alt[eq[:eq_Rref], endo[:Rref_t]] == ρ_smooth
    @test C_alt[eq[:eq_Rref]]                 == 0.
    @test Γ0_alt[eq[:eq_Rref], endo[:pgap_t]] == -φ_π * (1. - ρ_pgap) * (1. - ρ_smooth)
    @test Γ0_alt[eq[:eq_Rref], endo[:ygap_t]] == -φ_y * (1. - ρ_ygap) * (1. - ρ_smooth)
end

@testset "rw monetary policy rule: R_t = Rref_t + rw_t" begin
    # Unlike rw_zero_rate, R_t tracks the determinate reference rate (not a peg)
    @test Γ0_alt[eq[:eq_mp], endo[:R_t]]    == 1.
    @test Γ0_alt[eq[:eq_mp], endo[:Rref_t]] == -1.
    @test Γ0_alt[eq[:eq_mp], endo[:rw_t]]   == -1.
    @test C_alt[eq[:eq_mp]]                 == 0.
    @test all(Γ1_alt[eq[:eq_mp], :] .== 0.)
    @test all(Ψ_alt[eq[:eq_mp], :]  .== 0.)
    @test all(Π_alt[eq[:eq_mp], :]  .== 0.)
end

@testset "rw leaves the non-policy equilibrium conditions unchanged" begin
    @test Γ0_hist[eq_other, 1:n_st_hist] ≈ Γ0_alt[eq_other, 1:n_st_hist]
    @test Γ1_hist[eq_other, 1:n_st_hist] ≈ Γ1_alt[eq_other, 1:n_st_hist]
    @test Ψ_hist[eq_other, :]            ≈ Ψ_alt[eq_other, :]
    @test Π_hist[eq_other, :]            ≈ Π_alt[eq_other, :]
    @test !(Γ0_hist[old_eq_mp, 1:n_st_hist] ≈ Γ0_alt[old_eq_mp, 1:n_st_hist])
end

@testset "ait_Thalf / gdp_Thalf / ρ_rw settings feed through" begin
    mt = SmetsWouters(custom_settings = [Setting(:ait_Thalf, 20.),
                                         Setting(:gdp_Thalf, 5.),
                                         Setting(:ρ_rw, 0.5)])
    Γ0_t, Γ1_t, ~ = DSGE.rw_eqcond(mt)
    eqt, endot = mt.equilibrium_conditions, mt.endogenous_states
    @test Γ1_t[eqt[:eq_pgap], endot[:pgap_t]] == exp(log(0.5) / 20.)
    @test Γ1_t[eqt[:eq_ygap], endot[:ygap_t]] == exp(log(0.5) / 5.)
    @test Γ1_t[eqt[:eq_rw],   endot[:rw_t]]   == 0.5
end

@testset "rw_solve returns consistently sized, finite transition matrices" begin
    # rw is determinate (R_t tracks the reference rate), so gensys solves.
    TTT, RRR, CCC = DSGE.rw_solve(m)
    n = n_states_augmented(m)
    @test size(TTT) == (n, n)
    @test size(RRR) == (n, n_shocks_exogenous(m))
    @test length(CCC) == n
    @test all(isfinite, TTT)
    @test all(isfinite, RRR)
end

@testset "rw_forecast_init seeds pgap / ygap / rw / Rref" begin
    m <= Setting(:pgap_value, 2.5)
    m <= Setting(:ygap_value, 1.5)
    m <= Setting(:rw_value,   0.75)
    m <= Setting(:Rref_value, 0.1)

    n           = n_states_augmented(m)
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n)

    new_shocks, fs = DSGE.rw_forecast_init(m, shocks, copy(final_state))
    @test new_shocks == shocks
    @test length(fs) == n
    @test fs[endo[:pgap_t]] == -2.5
    @test fs[endo[:ygap_t]] == -1.5
    @test fs[endo[:rw_t]]   == -0.75
    @test fs[endo[:Rref_t]] == 0.1
    other = setdiff(1:n, [endo[:pgap_t], endo[:ygap_t], endo[:rw_t], endo[:Rref_t]])
    @test fs[other] == final_state[other]
end

@testset "rw_forecast_init falls back to R_t when :Rref_value unset" begin
    mf = SmetsWouters()
    DSGE.rw_eqcond(mf)
    mf <= Setting(:pgap_value, 0.)
    mf <= Setting(:ygap_value, 0.)
    mf <= Setting(:rw_value,   0.)
    # no :Rref_value
    n           = n_states_augmented(mf)
    shocks      = zeros(n_shocks_exogenous(mf), 1)
    final_state = collect(1.0:n)
    _, fs = DSGE.rw_forecast_init(mf, shocks, copy(final_state))
    @test fs[mf.endogenous_states[:Rref_t]] == final_state[mf.endogenous_states[:R_t]]
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    mb = SmetsWouters()
    mb <= Setting(:pgap_value, 2.5)
    mb <= Setting(:ygap_value, 1.5)
    mb <= Setting(:rw_value,   0.75)
    mb <= Setting(:Rref_value, 0.1)

    b_construct = @benchmark DSGE.rw()
    b_eqcond    = @benchmark DSGE.rw_eqcond($mb)
    b_solve     = @benchmark DSGE.rw_solve($mb)

    shocks      = zeros(n_shocks_exogenous(mb), 1)
    final_state = collect(1.0:n_states_augmented(mb))
    b_finit     = @benchmark DSGE.rw_forecast_init($mb, $shocks, copy($final_state))

    println("\n===== rw benchmark results =====")
    for (name, b) in [("rw",               b_construct),
                      ("rw_eqcond",        b_eqcond),
                      ("rw_solve",         b_solve),
                      ("rw_forecast_init", b_finit)]
        println(rpad(name, 20), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
