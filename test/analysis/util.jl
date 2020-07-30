using DSGE, Test, ModelConstructors

m = AnSchorfheide()
homedirpath = Sys.iswindows() ? joinpath(homedir(),".freddatarc") : joinpath(ENV["HOME"],".freddatarc")
if haskey(ENV, "FRED_API_KEY") || isfile(homedirpath)
    load_data(m)

    @testset "Test util functions" begin
        @test get_class(:histobs) == :obs
        @test get_class(:histpseudo) == :pseudo
        @test get_class(:histstates) == :states
        @test get_class(:histstdshock) == :stdshocks
        @test get_class(:histshocks) == :shocks
        @test_throws ErrorException get_class(:bad)

        @test get_product(:bddhistforecast4qobs) == :bddhistforecast4q
        @test get_product(:dettrendobs) == :dettrend
        @test_throws ErrorException get_product(:bad)

        @test DSGE.get_class_longname(:pseudo)==:pseudoobservable
        @test DSGE.get_class_longname(:obs)==:observable
        @test DSGE.get_class_longname(:states)==:state
        @test DSGE.get_class_longname(:shocks)==:shock
        @test DSGE.parse_transform(Symbol("DSGE.loggrowthtopct_annualized"))([1.0]) ==
            DSGE.loggrowthtopct_annualized([1.0])

        m <= Setting(:use_population_forecast, false)
        m <= Setting(:hpfilter_population, true)
        DSGE.load_population_growth(m)
        m <= Setting(:hpfilter_population, false)
        DSGE.load_population_growth(m)
        m <= Setting(:use_population_forecast, true)
        @test_throws ErrorException DSGE.load_population_growth(m)

        m <= Setting(:date_forecast_start, DSGE.quartertodate("2019-Q4"))
        m <= Setting(:date_presample_start, DSGE.quartertodate("1959-Q3"))
        @test DSGE.get_y0_index(m, :forecast) == 241
        @test DSGE.get_y0_index(m, :forecast4q) == 238
        @test DSGE.get_y0_index(m, :shockdec) == 2
        @test DSGE.get_y0_index(m, :hist) == 2
        @test DSGE.get_y0_index(m, :irf) == -1
        @test_throws ErrorException DSGE.get_y0_index(m, :bad)

        @test DSGE.check_consistent_order([1; 2; 3], [4; 5; 6])
        @test_throws AssertionError DSGE.check_consistent_order([1; 2; 3], [5; 4; 6])
        #would be nice to test:
        # resize population_forecast
        # get_population_series
        # get_mb_populuation_series

    end
else
    @warn "Skipping fred_data test because FRED_API_KEY not present"
end

@testset "Test prior_table works" begin
    DSGE.prior_table(m)

    DSGE.posterior_table(m, ones(16), ones(16, 2))

    DSGE.prior_posterior_moments_table(m, ones(16), ones(16, 2))
    DSGE.prior_posterior_table(m, ones(16))

end

@testset "Test finding density bands" begin
    @test find_density_bands(ones(100, 100), 0.95) == ones(2, 100)
    @test find_density_bands(ones(1, 100), 0.95) == ones(2, 100)
    @test find_density_bands(ones(100, 100), 0.95, minimize = true) == ones(2, 100)
    dens = find_density_bands(ones(100, 100), [.8, .9])
    @test propertynames(dens) == [Symbol("80.0% UB"), Symbol("80.0% LB"), Symbol("90.0% UB"), Symbol("90.0% LB")]
    @test Matrix(dens) == ones(100,4)
end
