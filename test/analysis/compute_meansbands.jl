using DSGE, HDF5, JLD
using Base.Test

path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:saveroot, tempdir())
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:forecast_uncertainty_override, Nullable(false))
m <= Setting(:use_population_forecast, true)
m <= Setting(:compute_shockdec_bands, true)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")

output_vars = add_requisite_output_vars([:histpseudo, :histobs,
                                         :hist4qpseudo, :hist4qobs,
                                         :forecastpseudo, :forecastobs,
                                         :forecast4qpseudo, :forecast4qobs,
                                         :shockdecpseudo, :shockdecobs,
                                         :irfpseudo, :irfobs])

# Read expected output
exp_modal_means, exp_modal_bands, exp_full_means, exp_full_bands =
    jldopen("$path/../reference/means_bands_out.jld", "r") do file
        read(file, "exp_modal_means"), read(file, "exp_modal_bands"),
        read(file, "exp_full_means"),  read(file, "exp_full_bands")
    end

# Modal
forecast_one(m, :mode, :none, output_vars, verbose = :none)
compute_meansbands(m, :mode, :none, output_vars; compute_shockdec_bands = true, verbose = :none)
meansbands_to_matrix(m, :mode, :none, output_vars; verbose = :none)

@testset "Check modal meansbands computation" begin
    for var in output_vars
        filename = get_forecast_filename(m, :mode, :none, Symbol("mb_matrix_", var),
                                         pathfcn = workpath, fileformat = :h5)
        @test @test_matrix_approx_eq exp_modal_means[var] h5read(filename, "means")
        @test @test_matrix_approx_eq exp_modal_bands[var] h5read(filename, "bands")
    end
end

# Full-distribution
@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
forecast_one(m, :full, :none, output_vars, verbose = :none)
compute_meansbands(m, :full, :none, output_vars; compute_shockdec_bands = true, verbose = :none)
meansbands_to_matrix(m, :full, :none, output_vars; verbose = :none)

@testset "Check full meansbands computation" begin
    for var in output_vars
        filename = get_forecast_filename(m, :full, :none, Symbol("mb_matrix_", var),
                                         pathfcn = workpath, fileformat = :h5)
        @test @test_matrix_approx_eq exp_full_means[var] h5read(filename, "means")
        @test @test_matrix_approx_eq exp_full_bands[var] h5read(filename, "bands")
    end
end

nothing
