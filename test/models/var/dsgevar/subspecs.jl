m1002 = Model1002("ss10"; custom_settings = [Setting(:add_laborshare_measurement, true),
                                             Setting(:add_NominalWageGrowth, true),
                                             Setting(:add_Epi_t_measurement, true)])
ans_dsge = AnSchorfheide()

@testset "DSGEVAR subspec ss0" begin
    m0 = DSGEVAR(m1002, collect(keys(m1002.exogenous_shocks)), "ss0")
    @test isempty(m0.observables)
    @test m0.lags == 0
    @test m0.λ == 0.
end

@testset "DSGEVAR subspecs ss1, ss10-ss13 (Model1002 observables)" begin
    m1  = DSGEVAR(m1002, collect(keys(m1002.exogenous_shocks)), "ss1")
    m10 = DSGEVAR(m1002, collect(keys(m1002.exogenous_shocks)), "ss10")
    m11 = DSGEVAR(m1002, collect(keys(m1002.exogenous_shocks)), "ss11")
    m12 = DSGEVAR(m1002, collect(keys(m1002.exogenous_shocks)), "ss12")
    m13 = DSGEVAR(m1002, collect(keys(m1002.exogenous_shocks)), "ss13")

    @test collect(keys(m1.observables)) == [:obs_hours, :obs_gdpdeflator]
    @test m1.lags == 4
    @test m1.λ == .5

    @test collect(keys(m10.observables)) == [:obs_hours, :obs_gdpdeflator,
                                              :laborshare_t, :NominalWageGrowth]
    @test m10.lags == 4
    @test m10.λ == .5

    @test collect(keys(m11.observables)) == [:obs_hours, :π_t,
                                              :laborshare_t, :NominalWageGrowth]
    @test m11.lags == 4
    @test m11.λ == .5

    @test collect(keys(m12.observables)) == [:obs_hours, :π_t,
                                              :laborshare_t, :NominalWageGrowth, :Epi_t]
    @test m12.lags == 4
    @test m12.λ == .5

    @test collect(keys(m13.observables)) == [:obs_spread, :obs_hours, :π_t,
                                              :laborshare_t, :NominalWageGrowth, :Epi_t]
    @test m13.lags == 4
    @test m13.λ == .5
end

@testset "DSGEVAR subspecs ss2-ss3 (AnSchorfheide observables)" begin
    m2 = DSGEVAR(ans_dsge, collect(keys(ans_dsge.exogenous_shocks)), "ss2")
    m3 = DSGEVAR(ans_dsge, collect(keys(ans_dsge.exogenous_shocks)), "ss3")

    @test collect(keys(m2.observables)) == [:obs_gdp, :obs_cpi]
    @test m2.lags == 4
    @test m2.λ == .5

    @test collect(keys(m3.observables)) == [:obs_gdp, :obs_cpi, :obs_nominalrate]
    @test m3.lags == 4
    @test m3.λ == .5
end

@testset "DSGEVAR undefined subspec" begin
    @test_throws ErrorException DSGEVAR(m1002, collect(keys(m1002.exogenous_shocks)), "ss999")
end

nothing
