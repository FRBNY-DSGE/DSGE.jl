using DSGE, Test, ModelConstructors, BenchmarkTools

m = AnSchorfheide()

pol = DSGE.default_policy()

@testset "default_policy AltPolicy object is wired up correctly" begin
    @test pol isa DSGE.AltPolicy
    @test pol.key == :default_policy
    # The default policy reuses the model's stock equilibrium conditions and solver
    @test pol.eqcond === eqcond
    @test pol.solve === solve
    # forecast_init defaults to identity (no adjustment to shocks / initial state)
    @test pol.forecast_init === identity
    # Constructor defaults
    @test pol.linestyle == :solid
end

@testset "default_policy eqcond reproduces the historical rule" begin
    Γ0_hist, Γ1_hist, C_hist, Ψ_hist, Π_hist = eqcond(m)
    Γ0_def,  Γ1_def,  C_def,  Ψ_def,  Π_def  = pol.eqcond(m)
    @test Γ0_def == Γ0_hist
    @test Γ1_def == Γ1_hist
    @test C_def  == C_hist
    @test Ψ_def  == Ψ_hist
    @test Π_def  == Π_hist
end

@testset "default_policy solve reproduces the historical rule" begin
    TTT_hist, RRR_hist, CCC_hist = solve(m)
    TTT_def,  RRR_def,  CCC_def  = pol.solve(m)
    @test TTT_def == TTT_hist
    @test RRR_def == RRR_hist
    @test CCC_def == CCC_hist
end

@testset "get_altpolicy(:default_policy) round-trips" begin
    pol2 = DSGE.get_altpolicy(:default_policy)
    @test pol2.key           == pol.key
    @test pol2.eqcond        === pol.eqcond
    @test pol2.solve         === pol.solve
    @test pol2.forecast_init === pol.forecast_init
end

################
# Benchmarking #
################
# Set this flag to true to run the default_policy benchmarks. Off by default so
# the test suite stays fast.
run_benchmarks = false

if run_benchmarks
    mb = AnSchorfheide()

    # Constructing the default_policy AltPolicy object
    b_construct = @benchmark DSGE.default_policy()

    # Looking it up through the get_altpolicy dispatcher
    b_lookup = @benchmark DSGE.get_altpolicy(:default_policy)

    # Evaluating the equilibrium conditions / solving via the wrapped functions
    b_eqcond = @benchmark $pol.eqcond($mb)
    b_solve  = @benchmark $pol.solve($mb)

    println("\n===== default_policy benchmark results =====")
    for (name, b) in [("default_policy", b_construct),
                      ("get_altpolicy",  b_lookup),
                      ("eqcond",         b_eqcond),
                      ("solve",          b_solve)]
        println(rpad(name, 16), " time: ", rpad(BenchmarkTools.prettytime(median(b).time), 12),
                "memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
