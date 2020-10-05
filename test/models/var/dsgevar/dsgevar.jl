dsge = Model1002("ss10")

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
    @test typeof(m.dsge) == typeof(dsge)
    @test spec(m) == "dsgevar_m1002"
    @test subspec(m) == "ss0"
    @test isempty(m.observables)
    @test !isempty(m.shocks)

    # Full construction w/shocks and model
    m = DSGE.DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss1")
    @test typeof(m.dsge) == typeof(dsge)
    @test spec(m) == "dsgevar_m1002"
    @test subspec(m) == "ss1"
    @test !isempty(m.observables)
    @test !isempty(m.shocks)

    # Check errors
    @test_throws ErrorException DSGE.DSGEVAR(dsge, [:not_a_shock], "ss0")
end

@testset "Helper functions for DSGE-VAR" begin

    # Check update
    as_dsge = AnSchorfheide()
    m = DSGE.DSGEVAR(as_dsge, collect(keys(as_dsge.exogenous_shocks)), "ss0")
    DSGE.update!(m; shocks = [:rm_sh])
    @test collect(keys(m.shocks)) == [:rm_sh] && isempty(m.observables) && m.lags == 0 && m.λ == 0.
    DSGE.update!(m; observables = [:obs_gdp])
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 0 && m.λ == 0.
    DSGE.update!(m; lags = 1)
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 1 && m.λ == 0.
    DSGE.update!(m; λ = 0.5)
    @test collect(keys(m.shocks)) == [:rm_sh] && collect(keys(m.observables)) == [:obs_gdp] && m.lags == 1 && m.λ == 0.5

    # Check errors
    @test_throws ErrorException DSGE.update!(m, dsge; shocks = [:not_a_shock])
    @test_throws ErrorException DSGE.update!(m, dsge; observables = [:not_an_observable])
    @test_throws AssertionError DSGE.update!(m, dsge; lags = -1)

    # Check access and size functions
    m = DSGE.DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss0")
    DSGE.update!(m; observables = [:obs_spread], lags = 1)
    @test collect(keys(DSGE.get_observables(m))) == [:obs_spread]
    @test n_observables(m) == 1
    @test collect(keys(DSGE.get_shocks(m))) == collect(keys(m.dsge.exogenous_shocks))
    @test DSGE.get_λ(m) == 0. == m.λ
    @test DSGE.n_shocks(m) == length(m.dsge.exogenous_shocks)
    @test DSGE.n_lags(m) == 1
    @test DSGE.get_lags(m) == 1
    @test typeof(DSGE.get_dsge(m)) == typeof(dsge)
end

fp = dirname(@__FILE__)
@testset "Pathing for DSGEVAR" begin
    ans_base = joinpath(fp, "save", "output_data", "m1002", "ss10", "dsgevar_m1002", "ss0", "estimate")
    ans_end = ["raw", "work", "tables", "figures", "log"]
    m = DSGE.DSGEVAR(dsge, Vector{Symbol}(undef, 0), "ss0")
    m <= Setting(:saveroot, joinpath(fp, "save"))
    m <= Setting(:dataroot, joinpath(fp, "save/input_data"))
    @test rawpath(m, "estimate") == joinpath(ans_base, ans_end[1])
    @test workpath(m, "estimate") == joinpath(ans_base, ans_end[2])
    @test tablespath(m, "estimate") == joinpath(ans_base, ans_end[3])
    @test figurespath(m, "estimate") == joinpath(ans_base, ans_end[4])
    @test logpath(m, "estimate") == joinpath(ans_base, ans_end[5])
end

@testset "DSGEVAR subspecs" begin
    dsge = Model1002("ss10"; custom_settings =
                  Dict{Symbol,Setting}(:add_laborshare_measurement =>
                                       Setting(:add_laborshare_measurement, true),
                                       :add_NominalWageGrowth =>
                                       Setting(:add_NominalWageGrowth, true)))
    m0  = DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss0")
    m1  = DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss1")
    m10 = DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss10")

    @test m0.λ  == 0.
    @test m1.λ  == 0.5
    @test m10.λ == 0.5
    @test m0.lags  == 0
    @test m1.lags  == 4
    @test m10.lags == 4
    @test isempty(m0.observables)
    @test collect(keys(m1.observables))  == [:obs_hours, :obs_gdpdeflator]
    @test collect(keys(m10.observables)) == [:obs_hours, :obs_gdpdeflator,
                                            :laborshare_t, :NominalWageGrowth]
    @test collect(values(m1.observables))  == collect(1:2)
    @test collect(values(m10.observables)) == collect(1:4)
end

nothing
