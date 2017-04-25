using DSGE, Base.Test, HDF5, JLD
include("../util.jl")

path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)
m <= Setting(:forecast_pseudoobservables, true)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")

output_vars = add_requisite_output_vars([:histpseudo, :histobs,
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
@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)


# Test read_forecast_output
for input_type in [:mode, :full]
    output_files = get_forecast_output_files(m, input_type, :none, output_vars)
    jldopen(output_files[:trendobs], "r") do file
        @test ndims(read_forecast_output(file, :obs, :trend, :obs_gdp)) == 2
    end
    jldopen(output_files[:forecastobs], "r") do file
        @test ndims(read_forecast_output(file, :obs, :forecast, :obs_gdp)) == 2
    end
    jldopen(output_files[:irfobs], "r") do file
        @test ndims(read_forecast_output(file, :obs, :irf, :obs_gdp, :rm_sh)) == 2
    end
end


nothing
