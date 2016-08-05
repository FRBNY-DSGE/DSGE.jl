using DSGE
using Base.Test

# Can we actually test? Require that the FRED_API_KEY ENV is populated.
if haskey(ENV, "FRED_API_KEY")
    # Specify exact vintage and dates
    custom_settings = Dict{Symbol, Setting}(
        :data_vintage         => Setting(:data_vintage, "151127"),
        :date_presample_start => Setting(:date_presample_start, quartertodate("2015-Q2")),
        :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
    m = Model990(custom_settings = custom_settings)
    m.testing=true

    df = load_data(m; try_disk=false, verbose=:none)

    # Silly comparisons
    delete!(df, :date)
    colmeans = colwise(mean, df)
    @test_approx_eq_eps colmeans[1][1] 0.4435 1e-4
    @test_approx_eq_eps colmeans[2][1] -51.2707 1e-4
    @test_approx_eq_eps colmeans[3][1] 0.1582 1e-4
    @test_approx_eq_eps colmeans[4][1] 0.4271 1e-4
    @test_approx_eq_eps colmeans[5][1] 0.3979 1e-4
    @test_approx_eq_eps colmeans[6][1] 0.0325 1e-4

    colstds  = colwise(std, df)
    @test_approx_eq_eps colstds[1][1] 0.3205 1e-4
    @test_approx_eq_eps colstds[2][1] 0.1015 1e-4
    @test_approx_eq_eps colstds[3][1] 0.3692 1e-4
    @test_approx_eq_eps colstds[4][1] 0.1367 1e-4
    @test_approx_eq_eps colstds[5][1] 0.0937 1e-4
    @test_approx_eq_eps colstds[6][1] 0.0000 1e-4
end
