path = dirname(@__FILE__)
fred = CSV.read("$path/../reference/fred_160812.csv")
custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :n_anticipated_shocks     => Setting(:n_anticipated_shocks, 6))

m = Model990(custom_settings = custom_settings, testing = true)
m <= Setting(:rate_expectations_source, :ois)
@testset "Check FRED data is properly loaded" begin
    @test @test_matrix_approx_eq Matrix(fred[:,2:end]) Matrix(load_fred_data(m, end_date = date_mainsample_end(m),
                                                                             verbose = :none)[:,2:end])
end
