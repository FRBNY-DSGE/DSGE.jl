using DSGE, Dates
path = dirname(@__FILE__)

mb_empty =  MeansBands()
@show mb_empty
@test isempty(mb_empty)
mb_full = load("$path/../reference/MeansBands.jld2", "mb")

@testset "Test that you can construct MeansBands objects and do stuff with them" begin
     # If one is empty, just return the non-empty
    @test cat(mb_empty, mb_full).means == mb_full.means
    @test cat(mb_full, mb_empty).means == mb_full.means
    @test get_class(mb_full) == :obs
    @test get_class(mb_empty) == :none
    @test get_product(mb_full) == :forecast
    @test get_product(mb_empty) == :none
    @test DSGE.get_cond_type(mb_full) == :none
    @test DSGE.get_cond_type(mb_empty) == :none
    @test DSGE.get_para(mb_full) == :full
    @test DSGE.get_para(mb_empty) == :none
    @test DSGE.n_vars_means(mb_full) == 3
    @test DSGE.n_vars_means(mb_empty) == 1
    @test DSGE.get_vars_means(mb_empty) == [:none]
    @test DSGE.get_vars_means(mb_full) == [:obs_gdp, :obs_cpi, :obs_nominalrate]
    @test DSGE.n_periods_means(mb_full) == 60
    @test DSGE.n_periods_means(mb_empty) == 1
    @test DSGE.startdate_means(mb_full) == quartertodate("2015-Q4")
    @test DSGE.startdate_means(mb_empty) == Dates.Date(0000, 1, 1)
    @test DSGE.enddate_means(mb_full) == quartertodate("2030-Q3")
    @test DSGE.enddate_means(mb_empty) == Dates.Date(0000, 1, 1)

    @test DSGE.n_vars_bands(mb_full) == 3
    @test DSGE.n_vars_bands(mb_empty) == 1
    @test DSGE.n_periods_bands(mb_full) == 60
    @test DSGE.n_periods_bands(mb_empty) == 1
    @test DSGE.startdate_bands(mb_full) == quartertodate("2015-Q4")
    @test DSGE.startdate_bands(mb_empty) == Dates.Date(0000, 1, 1)
    @test DSGE.enddate_bands(mb_full) == quartertodate("2030-Q3")
    @test DSGE.enddate_bands(mb_empty) == Dates.Date(0000, 1, 1)

    # Can only call on shockdec or irf MeansBands objects
    @test_throws AssertionError DSGE.get_shocks(mb_empty)
    @test_throws AssertionError DSGE.get_variables(mb_empty)
    @test_throws AssertionError DSGE.get_scenario_key(mb_empty)

    @test which_density_bands(mb_full) == ["90.0% LB", "80.0% LB", "70.0% LB", "60.0% LB", "50.0% LB", "50.0% UB", "60.0% UB", "70.0% UB", "80.0% UB", "90.0% UB"]
    @test which_density_bands(mb_full, uniquify = true) == ["50.0%", "60.0%", "70.0%", "80.0%", "90.0%"]
    @test which_density_bands(mb_empty) == String[]

    DSGE.prepare_meansbands_table_timeseries(mb_full, :obs_gdp)

    DSGE.parse_transform(Symbol("DSGE.loggrowthtopct_annualized"))([1.0]) == loggrowthtopct_annualized([1.0])

    # Ideally would also test prepare_means_table_shockdec...
end

@testset "Test meansbands_to_matrix works" begin
    meansbands_to_matrix(mb_full)
end
