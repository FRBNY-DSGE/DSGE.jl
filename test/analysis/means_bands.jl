using DSGE, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Initialize model object
custom_settings = Dict{Symbol, Setting}(
    :n_anticipated_shocks    => Setting(:n_anticipated_shocks, 6),
    :use_population_forecast => Setting(:use_population_forecast, true),
    :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :date_conditional_end    => Setting(:date_conditional_end, quartertodate("2015-Q4")),
    :forecast_kill_shocks    => Setting(:forecast_kill_shocks, true),
    :saveroot                => Setting(:saveroot, tempdir()))
m = Model990(custom_settings = custom_settings, testing = true)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference", "output_data", "m990", "ss2", "estimate", "raw"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "paramsmode_test.h5")
overrides[:full] = joinpath(estroot, "mhsave_test.h5")

output_vars = add_requisite_output_vars([:histpseudo,
                                         :forecastpseudo, :forecastobs,
                                         :shockdecpseudo, :shockdecobs,
                                         :irfpseudo, :irfobs])
@everywhere using DSGE

# Modal
@time forecast_one(m, :mode, :none, output_vars, verbose = :none)
@time means_bands_all(m, :mode, :none, output_vars; verbose = :none)
@time meansbands_matrix_all(m, :mode, :none, output_vars; verbose = :none)

# Full-distribution
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)
@time means_bands_all(m, :full, :none, output_vars; verbose = :none)
@time meansbands_matrix_all(m, :full, :none, output_vars; verbose = :none)

# TODO: Test means and bands against expected output


nothing
