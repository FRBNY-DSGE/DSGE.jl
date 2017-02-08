using DSGE, HDF5, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:forecast_kill_shocks, true)
m <= Setting(:use_population_forecast, true)
m <= Setting(:forecast_pseudoobservables, true)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")

output_vars = add_requisite_output_vars([:histpseudo,
                                         :forecastpseudo, :forecastobs,
                                         :shockdecpseudo, :shockdecobs,
                                         :irfpseudo, :irfobs])
@everywhere using DSGE

# Read expected output
exp_modal_means, exp_modal_bands, exp_full_means, exp_full_bands =
    jldopen("$path/../reference/means_bands_out.jld", "r") do file
        read(file, "exp_modal_means"), read(file, "exp_modal_bands"),
        read(file, "exp_full_means"),  read(file, "exp_full_bands")
    end

# Modal
@time forecast_one(m, :mode, :none, output_vars, verbose = :none)
@time means_bands_all(m, :mode, :none, output_vars; verbose = :none)
@time meansbands_matrix_all(m, :mode, :none, output_vars; verbose = :none)

mb_matrix_vars = map(x -> symbol("_matrix_$x"), output_vars)
files = get_meansbands_output_files(m, :mode, :none, mb_matrix_vars; fileformat = :h5)
for (var, mb_var) in zip(output_vars, mb_matrix_vars)
    filename = files[mb_var]
    @test_matrix_approx_eq exp_modal_means[var] h5read(filename, "means")
    @test_matrix_approx_eq exp_modal_bands[var] h5read(filename, "bands")
end

# Full-distribution
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)
@time means_bands_all(m, :full, :none, output_vars; verbose = :none)
@time meansbands_matrix_all(m, :full, :none, output_vars; verbose = :none)

files = get_meansbands_output_files(m, :full, :none, mb_matrix_vars; fileformat = :h5)
for (var, mb_var) in zip(output_vars, mb_matrix_vars)
    filename = files[mb_var]
    @test_matrix_approx_eq exp_full_means[var] h5read(filename, "means")
    @test_matrix_approx_eq exp_full_bands[var] h5read(filename, "bands")
end

nothing
