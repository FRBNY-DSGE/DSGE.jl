using DSGE, Base.Test, HDF5, JLD
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

# Check error handling for input_type = :subset
@test_throws ErrorException forecast_one(m, :subset, :none, output_vars,
                                subset_inds = 1:10, forecast_string = "",
                                verbose = :none)
@test_throws ErrorException forecast_one(m, :subset, :none, output_vars,
                                forecast_string = "test",
                                verbose = :none)

# Run modal forecasts
out = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
for cond_type in [:none, :semi, :full]
    @time forecast_one(m, :mode, cond_type, output_vars, verbose = :none)

    # Read output
    out[cond_type] = Dict{Symbol, Array{Float64}}()
    output_files = get_forecast_output_files(m, :mode, cond_type, output_vars)
    for var in keys(output_files)
        out[cond_type][var] = h5read(output_files[var], "arr")
    end
end

# Read expected output
exp_out = jldopen("$path/../reference/forecast_one_out.jld", "r") do file
    read(file, "exp_out")
end

# Test modal forecast outputs
specify_mode!(m, DSGE.get_forecast_input_file(m, :mode); verbose = :none)

for cond_type in [:none, :semi, :full]
    # Histories
    @test_matrix_approx_eq exp_out[cond_type][:histpseudo]     out[cond_type][:histpseudo]

    # Forecasts
    @test_matrix_approx_eq exp_out[cond_type][:forecastobs]    out[cond_type][:forecastobs]
    @test_matrix_approx_eq exp_out[cond_type][:forecastpseudo] out[cond_type][:forecastpseudo]

    # Shock decompositions, deterministic trends, trends
    @test_matrix_approx_eq exp_out[cond_type][:shockdecobs]    out[cond_type][:shockdecobs]
    @test_matrix_approx_eq exp_out[cond_type][:shockdecpseudo] out[cond_type][:shockdecpseudo]
    @test_matrix_approx_eq exp_out[cond_type][:dettrendobs]    out[cond_type][:dettrendobs]
    @test_matrix_approx_eq exp_out[cond_type][:dettrendpseudo] out[cond_type][:dettrendpseudo]
    @test_matrix_approx_eq exp_out[cond_type][:trendobs]       out[cond_type][:trendobs]
    @test_matrix_approx_eq exp_out[cond_type][:trendpseudo]    out[cond_type][:trendpseudo]

    # IRFs
    @test_matrix_approx_eq exp_out[cond_type][:irfobs]         out[cond_type][:irfobs]
    @test_matrix_approx_eq exp_out[cond_type][:irfpseudo]      out[cond_type][:irfpseudo]
end

# Test full-distribution blocking
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)


nothing
