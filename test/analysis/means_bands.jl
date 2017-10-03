using DSGE, HDF5, JLD

path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:forecast_uncertainty_override, Nullable(false))
m <= Setting(:use_population_forecast, true)
m <= Setting(:forecast_pseudoobservables, true)
m <= Setting(:compute_shockdec_bands, true)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")

output_vars = add_requisite_output_vars([:histpseudo, :histobs,
                                         :forecastpseudo, :forecastobs,
                                         :shockdecpseudo, :shockdecobs,
                                         :irfpseudo, :irfobs,
                                         :forecast4qobs, :bddforecast4qobs, :hist4qobs,
                                         :forecast4qpseudo, :bddforecast4qpseudo, :hist4qpseudo])

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

for var in output_vars
    filename = get_forecast_filename(m, :mode, :none, Symbol("mb_matrix_", var),
                                     pathfcn = workpath, fileformat = :h5)
    @test_matrix_approx_eq exp_modal_means[var] h5read(filename, "means")
    @test_matrix_approx_eq exp_modal_bands[var] h5read(filename, "bands")
end

# Full-distribution
@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)
@time means_bands_all(m, :full, :none, output_vars; verbose = :none)
@time meansbands_matrix_all(m, :full, :none, output_vars; verbose = :none)

for var in output_vars
    filename = get_forecast_filename(m, :full, :none, Symbol("mb_matrix_", var),
                                     pathfcn = workpath, fileformat = :h5)
    @test_matrix_approx_eq exp_full_means[var] h5read(filename, "means")
    @test_matrix_approx_eq exp_full_bands[var] h5read(filename, "bands")
end

nothing
