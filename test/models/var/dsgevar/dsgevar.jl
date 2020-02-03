using DSGE, Test, ModelConstructors, OrderedCollections

dsge = Model1002()

@testset "Constructing a DSGE-VAR" begin

    # Empty construction with a DSGE model
    m = DSGE.DSGEVAR(dsge, Vector{Symbol}(undef, 0), "ss0")
    @test typeof(m.dsge) == typeof(dsge)
    @test spec(m) == "dsgevar_m1002"
    @test subspec(m) == "ss0"
    @test isempty(m.observables)
    @test isempty(m.shocks)

    # Empty construction but w/some shocks
    m = DSGE.DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss0")

    # Full construction w/shocks and model


    # Check errors
    @test_throws ErrorException DSGE.DSGEVAR(dsge, [:not_a_shock], "ss0")
end

@testset "Helper functions for DSGE-VAR" begin

    # Check update
    as_dsge = AnSchorfheide()
    m = DSGE.DSGEVAR(as_dsge, collect(keys(as_dsge.exogenous_shocks)), "ss0")
    DSGE.update!(m; shocks = [:rm_sh])
    @test collect(keys(m.shocks)) == [:rm_sh] && isempty(m.observables) && m.lags == 0
    DSGE.update!(m; observables = [:obs_gdp])
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 0
    DSGE.update!(m; lags = 1)
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 1

    # Check errors
    @test_throws ErrorException DSGE.update!(m, dsge; shocks = [:not_a_shock])
    @test_throws ErrorException DSGE.update!(m, dsge; observables = [:not_an_observable])
    @test_throws AssertionError DSGE.update!(m, dsge; lags = -1)

    # Check access and size functions
    m = DSGE.DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss0")
    DSGE.update!(m; observables = [:obs_spread], lags = 1)
    @test DSGE.get_observables(m) == [:obs_spread]
    @test n_observables(m) == 1
    @test DSGE.get_shocks(m) == collect(keys(m.dsge.exogenous_shocks))
    @test DSGE.n_shocks(m) == length(m.dsge.exogenous_shocks)
    @test DSGE.n_lags(m) == 1
    @test typeof(DSGE.get_dsge(m)) == typeof(dsge)
end

nothing
