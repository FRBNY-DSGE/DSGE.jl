using DSGE
using StatsBase
using DataFrames
using HDF5, Base.Test

srand(1234)

customsettings = Dict{Symbol,Setting}(:date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = customsettings, testing = true)

testweights = [0.1, 0.3, 0.1, 0.01, 0.09, 0.15, 0.25]

test = multinomial_sampling2(m, testweights)
actual = [5.0;7.0;5.0;3.0;7.0;7.0;2.0]

@test_matrix_approx_eq test actual

nothing