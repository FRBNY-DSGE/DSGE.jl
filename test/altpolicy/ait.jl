using DSGE, Test, ModelConstructors, BenchmarkTools

m = AnSchorfheide()

# Record the historical-rule equilibrium conditions and indices BEFORE AIT
# augments the state space (ait_eqcond mutates the model by adding pgap).
Γ0_hist, Γ1_hist, ~ = eqcond(m)
n_eq_hist  = n_equilibrium_conditions(m)
n_st_hist  = n_states(m)
old_eq_mp  = m.equilibrium_conditions[:eq_mp]
eq_other   = setdiff(1:n_eq_hist, [old_eq_mp])   # non-policy rows
idx_R      = m.endogenous_states[:R_t]
idx_π      = m.endogenous_states[:π_t]

# Build the AIT equilibrium conditions. This adds :pgap_t / :eq_pgap.
Γ0_ait, Γ1_ait, C_ait, Ψ_ait, Π_ait = DSGE.ait_eqcond(m)

endo = m.endogenous_states
eq   = m.equilibrium_conditions

# Expected AIT coefficients (mirrors ait_replace_eq_entries)
ρ_ait = exp(log(0.5) / 10)
φ     = 0.25
ρ     = 0.0

@testset "AIT adds the pgap state and equilibrium condition" begin
    @test haskey(endo, :pgap_t)
    @test haskey(eq, :eq_pgap)
    @test endo[:pgap_t]  == n_st_hist + 1
    @test eq[:eq_pgap]   == n_eq_hist + 1
    # Augmented eqcond matrices are square at the new (larger) dimension
    @test size(Γ0_ait, 1) == size(Γ0_ait, 2) == n_states(m) == n_st_hist + 1
    @test size(Γ1_ait) == size(Γ0_ait)
end

@testset "AIT price-gap law of motion has the right coefficients" begin
    # pgap_t = π_t + ρ_ait * pgap_{t-1}
    @test Γ0_ait[eq[:eq_pgap], endo[:pgap_t]] == 1.
    @test Γ0_ait[eq[:eq_pgap], endo[:π_t]]    == -1.
    @test Γ1_ait[eq[:eq_pgap], endo[:pgap_t]] == ρ_ait
end

@testset "AIT monetary policy rule has the right coefficients" begin
    @test Γ0_ait[eq[:eq_mp], endo[:R_t]]    == 1.
    @test Γ0_ait[eq[:eq_mp], endo[:pgap_t]] == -φ * (1 / (1 - ρ_ait))
    @test Γ1_ait[eq[:eq_mp], endo[:R_t]]    == ρ
    @test C_ait[eq[:eq_mp]]                 == 0.
    # The old monetary policy rule's shock loadings are fully zeroed out
    @test Ψ_ait[eq[:eq_mp], :] == zeros(size(Ψ_ait, 2))
end

@testset "AIT leaves the non-policy equilibrium conditions unchanged" begin
    # Compare on the original (pre-augmentation) rows and columns
    @test Γ0_hist[eq_other, 1:n_st_hist] ≈ Γ0_ait[eq_other, 1:n_st_hist]
    @test Γ1_hist[eq_other, 1:n_st_hist] ≈ Γ1_ait[eq_other, 1:n_st_hist]
    # ...but the monetary policy row does change
    @test !(Γ0_hist[old_eq_mp, 1:n_st_hist] ≈ Γ0_ait[old_eq_mp, 1:n_st_hist])
end

@testset "ait_Thalf setting feeds through to ρ_ait" begin
    mt = AnSchorfheide(custom_settings = [Setting(:ait_Thalf, 20.)])
    Γ0_t, Γ1_t, ~ = DSGE.ait_eqcond(mt)
    @test Γ1_t[mt.equilibrium_conditions[:eq_pgap], mt.endogenous_states[:pgap_t]] ==
        exp(log(0.5) / 20.)
end

@testset "ait_forecast_init sets pgap to -pgap_value" begin
    pol         = DSGE.ait()
    n           = n_states_augmented(m)
    shocks      = zeros(n_shocks_exogenous(m), 1)
    final_state = collect(1.0:n)            # distinct values to detect mutation

    m <= Setting(:pgap_value, 0.)
    _, fs0 = pol.forecast_init(m, shocks, copy(final_state))
    @test fs0[endo[:pgap_t]] == 0.

    m <= Setting(:pgap_value, 2.5)
    _, fs1 = pol.forecast_init(m, shocks, copy(final_state))
    @test fs1[endo[:pgap_t]] == -2.5
    # Other states are left untouched
    other = setdiff(1:n, endo[:pgap_t])
    @test fs1[other] == final_state[other]
end

################
# Benchmarking #
################
# Set this flag to true to run the AIT benchmarks. Off by default so the test
# suite stays fast.
run_benchmarks = false

if run_benchmarks
    mb = Model990(custom_settings = [Setting(:add_altpolicy_pgap, true)])
    mb <= Setting(:pgap_type, :ait)
    mb <= Setting(:pgap_value, 0.)

    b_eqcond = @benchmark DSGE.ait_eqcond($mb)
    b_solve = @benchmark DSGE.ait_solve($mb)

    pol         = DSGE.ait()
    final_state = collect(1.0:n_states_augmented(mb))
    shocks      = zeros(n_shocks_exogenous(mb), 1)
    b_finit = @benchmark $pol.forecast_init($mb, $shocks, copy($final_state))

    println("\n===== AIT benchmark results =====")
    for (name, b) in [("ait_eqcond", b_eqcond), ("ait_solve", b_solve),
                      ("ait_forecast_init", b_finit)]
        println(rpad(name, 18), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
