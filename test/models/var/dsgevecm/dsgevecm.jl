dsge = Model1002("ss10")

@testset "Constructing a DSGE-VECM" begin

    # Empty construction with a DSGE model
    m = DSGE.DSGEVECM(dsge, Vector{Symbol}(undef, 0), "ss0")
    @test typeof(m.dsge) == typeof(dsge)
    @test spec(m) == "dsgevecm_m1002"
    @test subspec(m) == "ss0"
    @test isempty(m.observables)
    @test isempty(m.cointegrating)
    @test isempty(m.cointegrating_add)
    @test isempty(m.shocks)

    # Empty construction but w/some shocks
    m = DSGE.DSGEVECM(dsge, collect(keys(dsge.exogenous_shocks)), "ss0")
    @test typeof(m.dsge) == typeof(dsge)
    @test spec(m) == "dsgevecm_m1002"
    @test subspec(m) == "ss0"
    @test isempty(m.observables)
    @test isempty(m.cointegrating)
    @test isempty(m.cointegrating_add)
    @test !isempty(m.shocks)

    # Full construction w/shocks and model
    m = DSGE.DSGEVECM(dsge, collect(keys(dsge.exogenous_shocks)), "ss1")
    @test typeof(m.dsge) == typeof(dsge)
    @test spec(m) == "dsgevecm_m1002"
    @test subspec(m) == "ss1"
    @test !isempty(m.observables)
    @test isempty(m.cointegrating) # but cointegrating relationships still empty
    @test isempty(m.cointegrating_add)
    @test !isempty(m.shocks)

    # Check errors
    @test_throws ErrorException DSGE.DSGEVECM(dsge, [:not_a_shock], "ss0")
end

@testset "Helper functions for DSGE-VECM" begin

    # Check update
    as_dsge = AnSchorfheide()
    m = DSGE.DSGEVECM(as_dsge, collect(keys(as_dsge.exogenous_shocks)), "ss0")
    DSGE.update!(m; shocks = [:rm_sh])
    @test collect(keys(m.shocks)) == [:rm_sh] && isempty(m.observables) && m.lags == 0 && m.λ == 0. && isempty(m.cointegrating)
    DSGE.update!(m; observables = [:obs_gdp])
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 0 && m.λ == 0. && isempty(m.cointegrating)
    DSGE.update!(m; lags = 1)
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 1 && m.λ == 0. && isempty(m.cointegrating)
    DSGE.update!(m; λ = 0.5)
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 1 && m.λ == 0.5 && isempty(m.cointegrating)
    DSGE.update!(m; cointegrating = [:obs_gdp])
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 1 && m.λ == 0.5 && collect(keys(m.cointegrating)) == [:obs_gdp]
    DSGE.update!(m; cointegrating_add = [:obs_gdp])
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 1 && m.λ == 0.5 && collect(keys(m.cointegrating)) == [:obs_gdp] && collect(keys(m.cointegrating_add)) == [:obs_gdp]

    # Check errors
    @test_throws ErrorException DSGE.update!(m, dsge; shocks = [:not_a_shock])
    @test_throws ErrorException DSGE.update!(m, dsge; observables = [:not_an_observable])
    @test_throws ErrorException DSGE.update!(m, dsge; cointegrating = [:not_cointegrating])
    @test_throws AssertionError DSGE.update!(m, dsge; lags = -1)

    # Check access and size functions
    m = DSGE.DSGEVECM(dsge, collect(keys(dsge.exogenous_shocks)), "ss0")
    DSGE.update!(m; observables = [:obs_spread], lags = 1, cointegrating = [:obs_spread],
                 cointegrating_add = [:obs_spread])
    @test collect(keys(DSGE.get_observables(m))) == [:obs_spread]
    @test collect(keys(DSGE.get_cointegrating(m))) == [:obs_spread]
    @test collect(keys(DSGE.get_cointegrating(m))) == [:obs_spread]
    @test n_observables(m) == 1
    @test DSGE.n_cointegrating(m) == 1
    @test DSGE.n_cointegrating_add(m) == 1
    @test DSGE.get_observables(m)[:obs_spread] == 1
    @test DSGE.get_cointegrating(m)[:obs_spread] == 2
    @test DSGE.get_cointegrating_add(m)[:obs_spread] == 1
    @test collect(keys(DSGE.get_shocks(m))) == collect(keys(m.dsge.exogenous_shocks))
    @test DSGE.get_λ(m) == 0. == m.λ
    @test DSGE.n_shocks(m) == length(m.dsge.exogenous_shocks)
    @test DSGE.n_lags(m) == 1
    @test DSGE.get_lags(m) == 1
    @test typeof(DSGE.get_dsge(m)) == typeof(dsge)
end

nothing
