path = dirname(@__FILE__)

# Initialize model object
m = AnSchorfheide(testing = true)
m <= Setting(:cond_id, 0)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:use_population_forecast, true)

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")

skip_forecast_one_draw = false
if haskey(ENV, "FRED_API_KEY") || isfile(joinpath(homedir(),".freddatarc"))
    df = load_data(m)
    skip_forecast_one_draw = true
else
    @warn "Skipping forecast_one_draw tests because FRED_API_KEY not present"
end
# Make sure output_vars ignores the untransformed and 4Q things because they are
# computed in compute_meansbands
output_vars = add_requisite_output_vars([:histpseudo, :histobs, :histstdshocks,
                                         :histutpseudo, :histutobs,
                                         :hist4qpseudo, :hist4qobs,
                                         :forecaststates, :forecastpseudo, :forecastobs, :forecaststdshocks,
                                         :forecastutpseudo, :forecastutobs,
                                         :forecast4qpseudo, :forecast4qobs,
                                         :bddforecaststates, :bddforecastshocks, :bddforecastpseudo, :bddforecastobs,
                                         :shockdecpseudo, :shockdecobs,
                                         :trendstates, :trendobs, :trendpseudo,
                                         :dettrendstates, :dettrendobs, :dettrendpseudo,
                                         :irfstates, :irfpseudo, :irfobs])

# Check error handling for input_type = :subset
@testset "Ensure properly error handling for input_type = :subset" begin
    @test_throws ErrorException forecast_one(m, :subset, :none, output_vars,
                                    subset_inds = 1:10, forecast_string = "",
                                    verbose = :none)
    @test_throws ErrorException forecast_one(m, :subset, :none, output_vars,
                                    forecast_string = "test",
                                    verbose = :none)
end

# Run modal forecasts
out = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
for cond_type in [:none, :semi, :full]
    forecast_one(m, :mode, cond_type, output_vars, verbose = :none)
    forecast_one(m, :init, cond_type, output_vars, verbose = :none)

    # Read output
    out[cond_type] = Dict{Symbol, Array{Float64}}()
    output_files = get_forecast_output_files(m, :mode, cond_type, output_vars)
    for var in keys(output_files)
        out[cond_type][var] = load(output_files[var], "arr")
    end
end

# Read expected output
exp_out = JLD2.jldopen("$path/../reference/forecast_one_out.jld2", "r") do file
    read(file, "exp_out")
end

# Test modal forecast outputs
specify_mode!(m, DSGE.get_forecast_input_file(m, :mode); verbose = :none)

@testset "Test modal forecast for all major output_vars" begin
    for cond_type in [:none, :semi, :full]
        # Histories
        @test @test_matrix_approx_eq exp_out[cond_type][:histpseudo]     out[cond_type][:histpseudo]

        # Forecasts
        @test @test_matrix_approx_eq exp_out[cond_type][:forecastobs]    out[cond_type][:forecastobs]
        @test @test_matrix_approx_eq exp_out[cond_type][:forecastpseudo] out[cond_type][:forecastpseudo]

        # Shock decompositions, deterministic trends, trends
        @test @test_matrix_approx_eq exp_out[cond_type][:shockdecobs]    out[cond_type][:shockdecobs]
        @test @test_matrix_approx_eq exp_out[cond_type][:shockdecpseudo] out[cond_type][:shockdecpseudo]
        @test @test_matrix_approx_eq exp_out[cond_type][:dettrendobs]    out[cond_type][:dettrendobs]
        @test @test_matrix_approx_eq exp_out[cond_type][:dettrendpseudo] out[cond_type][:dettrendpseudo]
        @test @test_matrix_approx_eq exp_out[cond_type][:trendobs]       out[cond_type][:trendobs]
        @test @test_matrix_approx_eq exp_out[cond_type][:trendpseudo]    out[cond_type][:trendpseudo]

        # IRFs
        @test @test_matrix_approx_eq exp_out[cond_type][:irfobs]         out[cond_type][:irfobs]
        @test @test_matrix_approx_eq exp_out[cond_type][:irfpseudo]      out[cond_type][:irfpseudo]
    end
end

@testset "Test full-distribution forecasts run" begin
    m <= Setting(:forecast_block_size, 5)
    forecast_one(m, :prior, :none, output_vars, verbose = :none)
    for sampling_method in [:MH]
        m <= Setting(:sampling_method, sampling_method)
        for cond_type in [:none, :semi, :full]
            forecast_one(m, :full, cond_type, output_vars, verbose = :none)
            @test_throws ErrorException forecast_one(m, :subset, cond_type, output_vars, subset_inds = 1:10, verbose = :none)
            forecast_one(m, :subset, cond_type, output_vars, subset_inds = 1:10, forecast_string = "test", verbose = :none)
            forecast_one(m, :init_draw_shocks, cond_type, output_vars, verbose = :none)
            forecast_one(m, :mode_draw_shocks, cond_type, output_vars, verbose = :none)
        end
    end
end


# Test full-distribution blocking
@everywhere using DSGE
m <= Setting(:forecast_block_size, 5)
forecast_one(m, :full, :none, output_vars, verbose = :none)
# Test read_forecast_output
@testset "Test full-distribution blocking" begin
    for input_type in [:mode, :full]
        output_files = get_forecast_output_files(m, input_type, :none, output_vars)
        @test ndims(DSGE.read_forecast_series(output_files[:trendobs], :trend, m.observables[:obs_gdp])) == 2
        @test ndims(DSGE.read_forecast_series(output_files[:forecastobs], :forecast, m.observables[:obs_gdp])) == 2
        @test ndims(DSGE.read_forecast_series(output_files[:irfobs], m.observables[:obs_gdp], m.exogenous_shocks[:rm_sh])) == 2
    end
end

@testset "Test forecast_one_draw" begin
    if !skip_forecast_one_draw
        for input_type in [:mode, :full]
            params = if input_type == :mode
                load_draws(m, input_type)
            else
                load_draws(m, input_type)[1, :]
            end
            m <= Setting(:alternative_policy, AltPolicy(:historical, eqcond, solve))
            @test typeof(DSGE.forecast_one_draw(m, input_type, :none, output_vars, params, df)) == Dict{Symbol, Array{Float64}}

            # Test with alternative policy
            m <= Setting(:alternative_policy, AltPolicy(:taylor93, eqcond, solve))
            @test typeof(DSGE.forecast_one_draw(m, input_type, :none, output_vars, params, df)) == Dict{Symbol, Array{Float64}}
        end
    end
end

nothing
