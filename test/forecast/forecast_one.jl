using DSGE, Base.Test, JLD, DistributedArrays
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

output_vars = [:histpseudo, :forecastpseudo, :forecastobs, :bddforecastpseudo, :bddforecastobs,
               :shockdecpseudo, :shockdecobs, :dettrendpseudo, :dettrendobs, :trendpseudo, :trendobs,
               :irfpseudo, :irfobs]

# Check error handling for input_type = :subset
@test_throws ErrorException forecast_one(m, :subset, :none, output_vars,
                                subset_inds = 1:10, forecast_string = "",
                                verbose = :none)
@test_throws ErrorException forecast_one(m, :subset, :none, output_vars,
                                forecast_string = "test",
                                verbose = :none)

# Run modal forecasts
out = Dict{Symbol, Dict{Symbol, Any}}()
for cond_type in [:none, :semi, :full]
    @time dtemp = forecast_one(m, :mode, cond_type, output_vars, verbose = :none)

    # Convert to Dict of regular arrays
    temp = Dict{Symbol, Array}()
    for x in keys(dtemp)
        temp[x] = squeeze(convert(Array, dtemp[x]), 1)
    end
    out[cond_type] = temp
end

# Read expected output
exp_out = jldopen("$path/../reference/forecast_one_out.jld", "r") do file
    read(file, "exp_out")
end

# Test modal forecast outputs
specify_mode!(m, DSGE.get_input_file(m, :mode); verbose = :none)

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
m <= Setting(:forecast_blocking, false)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)

m <= Setting(:forecast_blocking, true)
@time forecast_one(m, :full, :none, output_vars, verbose = :none)

m <= Setting(:forecast_start_block, Nullable(2))
@test_throws ErrorException forecast_one(m, :full, :none, output_vars, verbose = :none)


nothing
