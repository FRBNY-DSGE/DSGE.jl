using DSGE, FileIO, JLD2, ModelConstructors, Test, Random, Dates, HDF5
path = dirname(@__FILE__)

generate_regime_switch_tests = false # Set to true if you want to regenerate the jld2 files for testing

if VERSION < v"1.5"
    ver = "111"
else
    ver = "150"
end

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

if haskey(ENV, "FRED_API_KEY") || isfile(joinpath(homedir(),".freddatarc"))
    df = load_data(m)
    skip_forecast_one_draw = false
else
    skip_forecast_one_draw = true
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
        m <= Setting(:regime_switching, true)
        for input_type in [:mode, :full]
            params = if input_type == :mode
                load_draws(m, input_type)
            else
                load_draws(m, input_type)[1, :]
            end
            setup_permanent_altpol!(m, AltPolicy(:historical, eqcond, solve))
            @test typeof(DSGE.forecast_one_draw(m, input_type, :none, output_vars, params, df;
                                                regime_switching = get_setting(m, :regime_switching),
                                                n_regimes = get_setting(m, :n_regimes))) == Dict{Symbol, Array{Float64}}

            # Test with alternative policy
            setup_permanent_altpol!(m, AltPolicy(:taylor93, eqcond, solve))
            @test typeof(DSGE.forecast_one_draw(m, input_type, :none, output_vars, params, df;
                                                regime_switching = get_setting(m, :regime_switching),
                                                n_regimes = get_setting(m, :n_regimes))) == Dict{Symbol, Array{Float64}}
        end
    end
end

## Now check for regime switching
# TODO: ADD TEST WHEN USING OVERRIDES W/REGIME-SWITCHING MODEL, also two more TODO below
custom_settings = Dict{Symbol, Setting}(
                                        :data_vintage             => Setting(:data_vintage, "160812"),
                                        :cond_vintage             => Setting(:cond_vintage, "160812"),
                                        :cond_id                  => Setting(:cond_id, 0),
                                        :use_population_forecast  => Setting(:use_population_forecast, true),
                                        :date_presample_start     => Setting(:date_presample_start, Date(1959, 9, 30)),
                                        :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
                                        :date_conditional_end     => Setting(:date_conditional_end, DSGE.quartertodate("2016-Q3")),
                                        :n_mon_anticipated_shocks => Setting(:n_mon_anticipated_shocks, 6))
:forecast_horizons => Setting(:forecast_horizons, 12)
:impulse_response_horizons => Setting(:impulse_response_horizons, 12)

# Now check that regime switching works for different types of possible regimes
regime_dates_dicts = [Dict{Int, Date}(1 => DSGE.quartertodate("1959-Q3"),
                                      2 => DSGE.quartertodate("2010-Q1"),
                                      3 => DSGE.quartertodate("2012-Q4")),
                      Dict{Int, Date}(1 => DSGE.quartertodate("1959-Q3"),
                                      2 => DSGE.quartertodate("1980-Q2"),
                                      3 => DSGE.quartertodate("2012-Q4")),
                      Dict{Int, Date}(1 => DSGE.quartertodate("1959-Q3"),
                                      2 => DSGE.quartertodate("1980-Q2"),
                                      3 => DSGE.quartertodate("2003-Q4")),
                      Dict{Int, Date}(1 => DSGE.quartertodate("1959-Q3"),
                                      2 => DSGE.quartertodate("2008-Q4"),
                                      3 => DSGE.quartertodate("2012-Q4")),
                      Dict{Int, Date}(1 => DSGE.quartertodate("1959-Q3"),
                                      2 => DSGE.quartertodate("2000-Q2"),
                                      3 => DSGE.quartertodate("2008-Q4"))]

# Read in the expected outputs for all the regimes except when the first one,
# which we use to test for shockdec and irf outputs. Those cases are handled separately
exp_out_dict = JLD2.jldopen("$path/../reference/forecast_one_out_rs1.jld2", "r") do file
    read(file, "exp_out_regime_switch_cases")
end

m = Model1002("ss10", custom_settings = custom_settings, testing = true)
m <= Setting(:rate_expectations_source, :ois)
dfs = Dict()
dfs[:none] = load("$path/../reference/regime_switch_data.jld2", "none")
dfs[:semi] = load("$path/../reference/regime_switch_data.jld2", "semi")
dfs[:full] = load("$path/../reference/regime_switch_data.jld2", "full")

if generate_regime_switch_tests
    exp_out_dict_new = exp_out_dict # We won't be testing, but we want to have the same structure as the existing one dict
end

@testset "Test modal and full distribution forecasts with regime switching for all major output_vars" begin
    # Loop over different times at which the regime switches, e.g. does the regime before or after ZLB
    m.settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
    m.test_settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
    m.settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
    m.test_settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
    for (k, v) in enumerate(regime_dates_dicts)

        # pseudo regime switching (no values have second/third regimes, but we are still saying regime_switching = true)
        m_rs1 = Model1002("ss10", testing = true, custom_settings = custom_settings)
        m_rs1 <= Setting(:rate_expectations_source, :ois)
        m_rs1.settings[:regime_switching] = Setting(:regime_switching, true)
        m_rs1.settings[:n_regimes] = Setting(:n_regimes, 3)
        m_rs1.settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
        m_rs1.settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
        m_rs1.settings[:regime_dates] = Setting(:regime_dates, v)
        m_rs1.settings[:time_varying_trends] = Setting(:time_varying_trends, false)
        m_rs1.test_settings[:regime_switching] = Setting(:regime_switching, true)
        m_rs1.test_settings[:n_regimes] = Setting(:n_regimes, 3)
        m_rs1.test_settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
        m_rs1.test_settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
        m_rs1.test_settings[:regime_dates] = Setting(:regime_dates, v)
        m_rs1.test_settings[:time_varying_trends] = Setting(:time_varying_trends, false)
        m_rs1 = DSGE.setup_regime_switching_inds!(m_rs1)

        # pseudo regime switching (identical values for standard deviations)
        m_rs2 = Model1002("ss51v", testing = true, custom_settings = custom_settings)
        m_rs2 <= Setting(:rate_expectations_source, :ois)
        m_rs2.settings[:regime_switching] = Setting(:regime_switching, true)
        m_rs2.settings[:n_regimes] = Setting(:n_regimes, 3)
        m_rs2.settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
        m_rs2.settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
        m_rs2.settings[:regime_dates] = Setting(:regime_dates, v)
        m_rs2.settings[:time_varying_trends] = Setting(:time_varying_trends, false)
        m_rs2.test_settings[:regime_switching] = Setting(:regime_switching, true)
        m_rs2.test_settings[:n_regimes] = Setting(:n_regimes, 3)
        m_rs2.test_settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
        m_rs2.test_settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
        m_rs2.test_settings[:time_varying_trends] = Setting(:time_varying_trends, false)
        m_rs2.test_settings[:regime_dates] = Setting(:regime_dates, v)
        m_rs2 = DSGE.setup_regime_switching_inds!(m_rs2)

        # non-trivial regime switching: eqcond matrices different and QQ matrix is different
        # - Does not test for switching in ZZ or ZZ_pseudo
        m_rs3 = Model1002("ss51v", testing = true, custom_settings = custom_settings)
        m_rs3 <= Setting(:rate_expectations_source, :ois)
        m_rs3.settings[:regime_switching] = Setting(:regime_switching, true)
        m_rs3.settings[:n_regimes] = Setting(:n_regimes, 3)
        m_rs3.settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
        m_rs3.settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
        m_rs3.settings[:regime_dates] = Setting(:regime_dates, v)
        m_rs3.settings[:time_varying_trends] = Setting(:time_varying_trends, false)
        m_rs3.test_settings[:regime_switching] = Setting(:regime_switching, true)
        m_rs3.test_settings[:n_regimes] = Setting(:n_regimes, 3)
        m_rs3.test_settings[:regime_switching_ndraws] = Setting(:regime_switching_ndraws, 4)
        m_rs3.test_settings[:forecast_ndraws] = Setting(:forecast_ndraws, 4)
        m_rs3.test_settings[:time_varying_trends] = Setting(:time_varying_trends, false)
        m_rs3.test_settings[:regime_dates] = Setting(:regime_dates, v)
        if haskey(m_rs3.settings, :model2para_regimes)
            delete!(m_rs3.settings, :model2para_regimes)
        end
        if haskey(m_rs3.test_settings, :model2para_regimes)
            delete!(m_rs3.test_settings, :model2para_regimes)
        end
        m_rs3 = DSGE.setup_regime_switching_inds!(m_rs3)

        # Need to set shocks for second and third regimes
        global prop = 1.   # Use the same value as the first regime
        global prop3 = .95 # Set the values to .95 of the first regime values
        for i in 1:3
            if i == 1
                oldprop3 = prop3
                prop3 = prop
            end
            ModelConstructors.set_regime_val!(m_rs2[:α], i, prop * m[:α].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_g], i, prop * m[:σ_g].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_b], i, prop * m[:σ_b].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_μ], i, prop * m[:σ_μ].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_ztil], i, prop * m[:σ_ztil].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_λ_f], i, prop * m[:σ_λ_f].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_λ_w], i, prop * m[:σ_λ_w].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_r_m], i, prop * m[:σ_r_m].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_σ_ω], i, prop * m[:σ_σ_ω].value)
            ModelConstructors.set_regime_val!(m_rs2[:σ_μ_e], i, prop * m[:σ_μ_e].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_γ], i, prop * m[:σ_γ].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_π_star], i, prop * m[:σ_π_star].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_lr], i, prop * m[:σ_lr].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_z_p], i, prop * m[:σ_z_p].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_tfp], i, prop * m[:σ_tfp].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_gdpdef], i, prop * m[:σ_gdpdef].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_corepce], i, prop * m[:σ_corepce].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_gdp], i, prop * m[:σ_gdp].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs2[:σ_gdi], i, prop * m[:σ_gdi].value, override_bounds = true)

            for j = 1:DSGE.n_mon_anticipated_shocks(m_rs2)
                ModelConstructors.set_regime_val!(m_rs2[Symbol("σ_r_m$(j)")], i, prop * m[Symbol("σ_r_m$(i)")])
            end

            ModelConstructors.set_regime_val!(m_rs3[:α], i, prop3 * m[:α].value) # Change α to change eqcond matrices
            ModelConstructors.set_regime_val!(m_rs3[:σ_g], i, prop3 * m[:σ_g].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_b], i, prop3 * m[:σ_b].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_μ], i, prop3 * m[:σ_μ].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_ztil], i, prop3 * m[:σ_ztil].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_λ_f], i, prop3 * m[:σ_λ_f].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_λ_w], i, prop3 * m[:σ_λ_w].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_r_m], i, prop3 * m[:σ_r_m].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_σ_ω], i, prop3 * m[:σ_σ_ω].value)
            ModelConstructors.set_regime_val!(m_rs3[:σ_μ_e], i, prop3 * m[:σ_μ_e].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_γ], i, prop3 * m[:σ_γ].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_π_star], i, prop3 * m[:σ_π_star].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_lr], i, prop3 * m[:σ_lr].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_z_p], i, prop3 * m[:σ_z_p].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_tfp], i, prop3 * m[:σ_tfp].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_gdpdef], i, prop3 * m[:σ_gdpdef].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_corepce], i, prop3 * m[:σ_corepce].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_gdp], i, prop3 * m[:σ_gdp].value, override_bounds = true)
            ModelConstructors.set_regime_val!(m_rs3[:σ_gdi], i, prop3 * m[:σ_gdi].value, override_bounds = true)

            for j = 1:DSGE.n_mon_anticipated_shocks(m_rs3)
                ModelConstructors.set_regime_val!(m_rs3[Symbol("σ_r_m$(j)")], i, prop3 * m[Symbol("σ_r_m$(i)")])
            end

            if i == 1
                global prop3 = oldprop3
            end
        end

        # Set up output vars, w/shockdecs and irfs tested only once b/c they take up a lot of space
        if k == 1
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
        else
            output_vars = add_requisite_output_vars([:histpseudo, :histobs, :histstdshocks,
                                                     :histutpseudo, :histutobs,
                                                     :hist4qpseudo, :hist4qobs,
                                                     :forecaststates, :forecastpseudo, :forecastobs, :forecaststdshocks,
                                                     :forecastutpseudo, :forecastutobs,
                                                     :forecast4qpseudo, :forecast4qobs,
                                                     :bddforecaststates, :bddforecastshocks, :bddforecastpseudo, :bddforecastobs])
        end

        # Save output in these dictionaries for different cases
        out = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
        out_rs1 = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
        out_rs2 = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
        out_rs3 = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
        params = map(x -> x.value, m.parameters)
        params_rs1 = map(x -> x.value, m_rs1.parameters)
        params_rs2 = map(x -> x.value, m_rs2.parameters)
        for para in m_rs3.parameters
            ModelConstructors.toggle_regime!(para, 3)
        end
        params_rs3 = map(x -> x.value, m_rs3.parameters)

        # Run modal forecasts w/different conditioning
        for cond_type in [:none, :semi, :full]
            setup_regime_switching_inds!(m_rs1; cond_type = cond_type)
            setup_regime_switching_inds!(m_rs2; cond_type = cond_type)
            setup_regime_switching_inds!(m_rs3; cond_type = cond_type)
            forecast_one(m, :mode, cond_type, output_vars, verbose = :none, params = params, df = dfs[cond_type])
            forecast_one(m, :init, cond_type, output_vars, verbose = :none, params = params, df = dfs[cond_type])
            forecast_one(m_rs1, :mode, cond_type, output_vars, verbose = :none, params = params_rs1, df = dfs[cond_type])
            forecast_one(m_rs1, :init, cond_type, output_vars, verbose = :none, params = params_rs1, df = dfs[cond_type])
            forecast_one(m_rs2, :mode, cond_type, output_vars, verbose = :none, params = params_rs2, df = dfs[cond_type])
            forecast_one(m_rs2, :init, cond_type, output_vars, verbose = :none, params = params_rs2, df = dfs[cond_type])
            forecast_one(m_rs3, :mode, cond_type, output_vars, verbose = :none, params = params_rs3, df = dfs[cond_type])
            forecast_one(m_rs3, :init, cond_type, output_vars, verbose = :none, params = params_rs3, df = dfs[cond_type])

            # Read output
            out[cond_type] = Dict{Symbol, Array{Float64}}()
            out_rs1[cond_type] = Dict{Symbol, Array{Float64}}()
            out_rs2[cond_type] = Dict{Symbol, Array{Float64}}()
            out_rs3[cond_type] = Dict{Symbol, Array{Float64}}()
            output_files = get_forecast_output_files(m, :mode, cond_type, output_vars)
            output_files_rs1 = get_forecast_output_files(m_rs1, :mode, cond_type, output_vars)
            output_files_rs2 = get_forecast_output_files(m_rs2, :mode, cond_type, output_vars)
            output_files_rs3 = get_forecast_output_files(m_rs3, :mode, cond_type, output_vars)
            for var in keys(output_files)
                if !(var in [:shockdecpseudo, :shockdecobs, :irfstates, :irfobs, :irfpseudo]) ||
                    (k == 1)
                    out[cond_type][var] = load(output_files[var], "arr")
                    out_rs1[cond_type][var] = load(output_files_rs1[var], "arr")
                    out_rs2[cond_type][var] = load(output_files_rs2[var], "arr")
                    out_rs3[cond_type][var] = load(output_files_rs3[var], "arr")
                end
            end
        end

        # Check parameters actually switched and can switch back
        for para in m_rs3.parameters
            ModelConstructors.toggle_regime!(para, 3)
        end
        @test m_rs3[:α].value != m[:α].value # check the regimes do not match
        for para in m_rs3.parameters
            ModelConstructors.toggle_regime!(para, 1)
        end
        @test m_rs3[:α].value == m[:α].value # check the regimes match after toggling

        # Read expected output
        if k == 1
            exp_out, exp_out_true = JLD2.jldopen("$path/../reference/forecast_one_out_rs2_version=" * ver * ".jld2", "r") do file
                read(file, "exp_out_regime_switch"), read(file, "exp_out_true_regime_switch")
            end
        else
            exp_out      = exp_out_dict[k][:out]
            exp_out_true = exp_out_dict[k][:out_rs3]
        end

        if generate_regime_switch_tests
            if k == 1
                exp_out_regime_switch_new = out
                exp_out_true_regime_switch_new = out_rs3
            else
                exp_out_dict_new[k][:out] = out
                exp_out_dict_new[k][:out_rs3] = out_rs3
            end
        else
            for cond_type in [:none]#, :semi, :full]
                # Histories
                @test @test_matrix_approx_eq exp_out[cond_type][:histpseudo]          out[cond_type][:histpseudo]
                @test @test_matrix_approx_eq exp_out[cond_type][:histpseudo]          out_rs1[cond_type][:histpseudo]
                @test @test_matrix_approx_eq exp_out[cond_type][:histpseudo]          out_rs2[cond_type][:histpseudo]
                @test @test_matrix_approx_eq exp_out_true[cond_type][:histpseudo]     out_rs3[cond_type][:histpseudo]
                @test !(exp_out[cond_type][:histpseudo] ≈                             out_rs3[cond_type][:histpseudo])

                # Forecasts
                @test @test_matrix_approx_eq exp_out[cond_type][:forecastobs]             out[cond_type][:forecastobs]
                @test @test_matrix_approx_eq exp_out[cond_type][:forecastpseudo]          out[cond_type][:forecastpseudo]

                @test @test_matrix_approx_eq exp_out[cond_type][:forecastobs]             out_rs1[cond_type][:forecastobs]
                @test @test_matrix_approx_eq exp_out[cond_type][:forecastobs]             out_rs2[cond_type][:forecastobs]
                @test @test_matrix_approx_eq exp_out_true[cond_type][:forecastobs]        out_rs3[cond_type][:forecastobs]
                @test @test_matrix_approx_eq exp_out[cond_type][:forecastpseudo]          out_rs1[cond_type][:forecastpseudo]
                @test @test_matrix_approx_eq exp_out[cond_type][:forecastpseudo]          out_rs2[cond_type][:forecastpseudo]
                @test @test_matrix_approx_eq exp_out_true[cond_type][:forecastpseudo]     out_rs3[cond_type][:forecastpseudo]
                @test !(exp_out[cond_type][:forecastobs] ≈                                out_rs3[cond_type][:forecastobs])
                @test !(exp_out[cond_type][:forecastpseudo] ≈                             out_rs3[cond_type][:forecastpseudo])

                # B/c memory limits, only testing shock decs w/first test case
                if k == 1
                    # Shock decompositions, deterministic trends, trends
                    @test @test_matrix_approx_eq exp_out[cond_type][:shockdecobs]    out[cond_type][:shockdecobs]
                    @test @test_matrix_approx_eq exp_out[cond_type][:shockdecpseudo] out[cond_type][:shockdecpseudo]
                    @test @test_matrix_approx_eq exp_out[cond_type][:dettrendobs]    out[cond_type][:dettrendobs]
                    @test @test_matrix_approx_eq exp_out[cond_type][:dettrendpseudo] out[cond_type][:dettrendpseudo]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendobs]       out[cond_type][:trendobs]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendpseudo]    out[cond_type][:trendpseudo]

                    @test @test_matrix_approx_eq exp_out[cond_type][:shockdecobs]    out_rs1[cond_type][:shockdecobs]
                    @test @test_matrix_approx_eq exp_out[cond_type][:shockdecpseudo] out_rs1[cond_type][:shockdecpseudo]
                    @test @test_matrix_approx_eq exp_out[cond_type][:dettrendobs]    out_rs1[cond_type][:dettrendobs]
                    @test @test_matrix_approx_eq exp_out[cond_type][:dettrendpseudo] out_rs1[cond_type][:dettrendpseudo]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendobs]       out_rs1[cond_type][:trendobs][:, 1]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendpseudo]    out_rs1[cond_type][:trendpseudo][:, 1]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendobs]       out_rs1[cond_type][:trendobs][:, 2]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendpseudo]    out_rs1[cond_type][:trendpseudo][:, 2]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendobs]       out_rs1[cond_type][:trendobs][:, 3]
                    @test @test_matrix_approx_eq exp_out[cond_type][:trendpseudo]    out_rs1[cond_type][:trendpseudo][:, 3]

                    @test @test_matrix_approx_eq exp_out_true[cond_type][:shockdecobs]    out_rs3[cond_type][:shockdecobs]
                    @test @test_matrix_approx_eq exp_out_true[cond_type][:shockdecpseudo] out_rs3[cond_type][:shockdecpseudo]
                    @test @test_matrix_approx_eq exp_out_true[cond_type][:dettrendobs]    out_rs3[cond_type][:dettrendobs]
                    @test @test_matrix_approx_eq exp_out_true[cond_type][:dettrendpseudo] out_rs3[cond_type][:dettrendpseudo]
                    @test @test_matrix_approx_eq exp_out_true[cond_type][:trendobs]       out_rs3[cond_type][:trendobs]
                    @test @test_matrix_approx_eq exp_out_true[cond_type][:trendpseudo]    out_rs3[cond_type][:trendpseudo]

                    # IRFs
                    @test @test_matrix_approx_eq exp_out[cond_type][:irfobs]              out[cond_type][:irfobs]
                    @test @test_matrix_approx_eq exp_out[cond_type][:irfpseudo]           out[cond_type][:irfpseudo]
                    @test @test_matrix_approx_eq exp_out[cond_type][:irfobs]              out_rs1[cond_type][:irfobs]
                    @test @test_matrix_approx_eq exp_out[cond_type][:irfpseudo]           out_rs1[cond_type][:irfpseudo]
                    @test @test_matrix_approx_eq exp_out_true[cond_type][:irfobs]         out_rs3[cond_type][:irfobs]
                    @test @test_matrix_approx_eq exp_out_true[cond_type][:irfpseudo]      out_rs3[cond_type][:irfpseudo]
                end
            end
        end

        if k == 1 # only testing full distribution with the first case of regime switching

            # Construct fake matrix of parameters draws
            params_rs1 = repeat(params_rs1', get_setting(m_rs1, :regime_switching_ndraws))
            params_rs2 = repeat(params_rs2', get_setting(m_rs2, :regime_switching_ndraws))
            params_rs3 = repeat(params_rs3', get_setting(m_rs3, :regime_switching_ndraws))
            m <= Setting(:forecast_block_size, 2)
            m_rs1 <= Setting(:forecast_block_size, 2)
            m_rs2 <= Setting(:forecast_block_size, 2)
            m_rs3 <= Setting(:forecast_block_size, 2)
            m <= Setting(:forecast_jstep, 1)
            m_rs1 <= Setting(:forecast_jstep, 1)
            m_rs2 <= Setting(:forecast_jstep, 1)
            m_rs3 <= Setting(:forecast_jstep, 1)

            # @testset "Test full-distribution forecasts run" begin
            # TODO: UNCOMMENT THE FOLLOWING LINES WHEN LOADING DRAWS FROM A SAVED ESTIMATION HAS BEEN IMPLEMENTED
            # for sampling_method in [:MH]
            #     m <= Setting(:sampling_method, sampling_method)
            #     for cond_type in [:none, :semi, :full]
            #         forecast_one(m, :full, cond_type, output_vars, verbose = :none)
            #         @test_throws ErrorException forecast_one(m, :subset, cond_type, output_vars, subset_inds = 1:10, verbose = :none)
            #         forecast_one(m, :subset, cond_type, output_vars, subset_inds = 1:10, forecast_string = "test", verbose = :none)
            #         forecast_one(m, :init_draw_shocks, cond_type, output_vars, verbose = :none)
            #         forecast_one(m, :mode_draw_shocks, cond_type, output_vars, verbose = :none)
            #     end
            # end
            # end

            out_rs1 = Dict{Symbol, Dict{Symbol, Dict{Symbol, Array{Float64}}}}()
            out_rs2 = Dict{Symbol, Dict{Symbol, Dict{Symbol, Array{Float64}}}}()
            out_rs3 = Dict{Symbol, Dict{Symbol, Dict{Symbol, Array{Float64}}}}()

            exp_out = JLD2.jldopen("$path/../reference/forecast_one_out_rs2_version=" * ver * ".jld2", "r") do file
                read(file, "exp_out_regime_switch_full")
            end
            exp_out_true = JLD2.jldopen("$path/../reference/forecast_one_out_rs3_version=" * ver * ".jld2", "r") do file
                read(file, "exp_out_true_regime_switch_full")
            end

            for cond_type in [:none, :semi, :full]
                # out[cond_type] = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
                out_rs1[cond_type] = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
                out_rs2[cond_type] = Dict{Symbol, Dict{Symbol, Array{Float64}}}()
                out_rs3[cond_type] = Dict{Symbol, Dict{Symbol, Array{Float64}}}()

                # Forecast and read output from forecast
                for fcast_type in [:full, :subset, :init_draw_shocks, :mode_draw_shocks, :prior]
                    if fcast_type == :subset
                        for (model, para_rs) in zip([m_rs1, m_rs2, m_rs3], [params_rs1, params_rs2, params_rs3])
                            Random.seed!(1793)
                            forecast_one(model, :subset, cond_type, output_vars, subset_inds = 1:2, forecast_string = "test",
                                         verbose = :none, params = para_rs, df = dfs[cond_type])
                            @test_throws ErrorException forecast_one(model, :subset, cond_type, output_vars, subset_inds = 1:2,
                                                                     verbose = :none, params = para_rs, df = dfs[cond_type])
                        end
                    elseif cond_type == :none && fcast_type == :prior
                        # TODO: IMPLEMENT DRAWING FROM PRIOR FOR FORECAST WITH REGIME SWITCHING
                        # forecast_one(m_rs1, :prior, :none, output_vars, verbose = :none)
                        continue
                    elseif fcast_type == :prior
                        continue
                    elseif fcast_type == :mode_draw_shocks
                        for (model, para_rs) in zip([m_rs1, m_rs2, m_rs3], [params_rs1, params_rs2, params_rs3])
                            Random.seed!(1793)
                            forecast_one(model, fcast_type, cond_type, output_vars, verbose = :none, params = para_rs[1, :], df = dfs[cond_type])
                        end
                    else
                        for (model, para_rs) in zip([m_rs1, m_rs2, m_rs3], [params_rs1, params_rs2, params_rs3])
                            Random.seed!(1793)
                            forecast_one(model, fcast_type, cond_type, output_vars, verbose = :none, params = para_rs, df = dfs[cond_type])
                        end
                    end
                    out_rs1[cond_type][fcast_type] = Dict{Symbol, Array{Float64}}()
                    out_rs2[cond_type][fcast_type] = Dict{Symbol, Array{Float64}}()
                    out_rs3[cond_type][fcast_type] = Dict{Symbol, Array{Float64}}()
                    if fcast_type == :subset
                        output_files_rs1 = get_forecast_output_files(m_rs1, fcast_type, cond_type, output_vars; forecast_string = "test")
                        output_files_rs2 = get_forecast_output_files(m_rs2, fcast_type, cond_type, output_vars; forecast_string = "test")
                        output_files_rs3 = get_forecast_output_files(m_rs3, fcast_type, cond_type, output_vars; forecast_string = "test")
                    else
                        output_files_rs1 = get_forecast_output_files(m_rs1, fcast_type, cond_type, output_vars)
                        output_files_rs2 = get_forecast_output_files(m_rs2, fcast_type, cond_type, output_vars)
                        output_files_rs3 = get_forecast_output_files(m_rs3, fcast_type, cond_type, output_vars)
                    end
                    for var in keys(output_files_rs1)
                        if fcast_type != :prior || (cond_type == :none)
                            out_rs1[cond_type][fcast_type][var] = load(output_files_rs1[var], "arr")
                        end
                        if fcast_type != :prior
                            out_rs2[cond_type][fcast_type][var] = load(output_files_rs2[var], "arr")
                            out_rs3[cond_type][fcast_type][var] = load(output_files_rs3[var], "arr")
                        end
                    end
                end
            end

            if generate_regime_switch_tests
                exp_out_regime_switch_full_new = out_rs1
                exp_out_true_regime_switch_full_new = out_rs3

                JLD2.jldopen("$path/../reference/forecast_one_out_rs2_version=" * ver * ".jld2", true, true, true, IOStream) do file
                    # Now for the new additions
                    write(file, "exp_out_regime_switch", exp_out_regime_switch_new)
                    write(file, "exp_out_true_regime_switch", exp_out_true_regime_switch_new)
                    write(file, "exp_out_regime_switch_full", exp_out_regime_switch_full_new)
                end

                JLD2.jldopen("$path/../reference/forecast_one_out_rs3_version=" * ver * ".jld2", true, true, true, IOStream) do file
                    write(file, "exp_out_true_regime_switch_full", exp_out_true_regime_switch_full_new)
                end
            else
                for cond_type in [:none, :semi, :full]
                    for fcast_type in [:full, :subset, :init_draw_shocks, :mode_draw_shocks] # TODO: ADD PRIOR
                        # Histories
                        @test maximum(abs.(exp_out[cond_type][fcast_type][:histpseudo] - out_rs1[cond_type][fcast_type][:histpseudo])) < 1e-5
                        @test maximum(abs.(exp_out[cond_type][fcast_type][:histpseudo] - out_rs2[cond_type][fcast_type][:histpseudo])) < 1e-5
                        @test maximum(abs.(exp_out_true[cond_type][fcast_type][:histpseudo] - out_rs3[cond_type][fcast_type][:histpseudo])) < 1e-5
                        @test !(exp_out[cond_type][fcast_type][:histpseudo] ≈                             out_rs3[cond_type][fcast_type][:histpseudo])
                        # Forecasts
                        @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:forecastobs]             out_rs1[cond_type][fcast_type][:forecastobs]
                        @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:forecastobs]             out_rs2[cond_type][fcast_type][:forecastobs]
                        @test @test_matrix_approx_eq exp_out_true[cond_type][fcast_type][:forecastobs][:, vcat(1:9, 13), 1]        out_rs3[cond_type][fcast_type][:forecastobs][:, vcat(1:9, 13), 1]
                        @test @test_matrix_approx_eq exp_out_true[cond_type][fcast_type][:forecastobs][:, :, 2]        out_rs3[cond_type][fcast_type][:forecastobs][:, :, 2]
                        @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:forecastpseudo]          out_rs1[cond_type][fcast_type][:forecastpseudo]
                        @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:forecastpseudo]          out_rs2[cond_type][fcast_type][:forecastpseudo]
                        @test @test_matrix_approx_eq exp_out_true[cond_type][fcast_type][:forecastpseudo]     out_rs3[cond_type][fcast_type][:forecastpseudo]
                        @test !(exp_out[cond_type][fcast_type][:forecastobs] ≈                                out_rs3[cond_type][fcast_type][:forecastobs])
                        @test !(exp_out[cond_type][fcast_type][:forecastpseudo] ≈                             out_rs3[cond_type][fcast_type][:forecastpseudo])

                        # We do not compare distributional shocks decs or IRFs b/c memory limits, although we do check that the methods run
                        if k == 1
                            # # Shock decompositions, deterministic trends, trends
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:shockdecobs]    out[cond_type][fcast_type][:shockdecobs]
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:shockdecpseudo] out[cond_type][fcast_type][:shockdecpseudo]
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:dettrendobs]    out[cond_type][fcast_type][:dettrendobs]
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:dettrendpseudo] out[cond_type][fcast_type][:dettrendpseudo]
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:trendobs]       out[cond_type][fcast_type][:trendobs]
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:trendpseudo]    out[cond_type][fcast_type][:trendpseudo]

                            # # IRFs
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:irfobs]         out[cond_type][fcast_type][:irfobs]
                            # @test @test_matrix_approx_eq exp_out[cond_type][fcast_type][:irfpseudo]      out[cond_type][fcast_type][:irfpseudo]
                        end
                    end
                end
            end
        end
    end

    if generate_regime_switch_tests
        JLD2.jldopen("$path/../reference/forecast_one_out_rs1_version=" * ver * ".jld2", true, true, true, IOStream) do file

            # Now for the new additions
            write(file, "exp_out_regime_switch_cases", exp_out_dict_new)
        end
    end

end


nothing
