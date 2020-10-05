using DSGE, Test, ModelConstructors, Plots

if haskey(ENV, "FRED_API_KEY")

    GR.inline("pdf")
    GR.inline("png")

    m = Model1002("ss10")

    usual_model_settings!(m, "191010")
    m <= Setting(:use_population_forecast, false)
    estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
    overrides = forecast_input_file_overrides(m)
    overrides[:mode] = joinpath(estroot, "optimize_1002.h5")
    overrides[:full] = joinpath(estroot, "metropolis_hastings_1002.h5")

    @testset "Ensure packet drivers run without deprecation" begin
        usual_model_forecast(m, :mode, :none, [:histobs, :histpseudo, :histstates, :forecastobs, :forecastpseudo, :forecaststates, :forecast4qobs, :forecast4qpseudo], mb_matrix = true, check_empty_columns = false)
        m <= Setting(:forecast_jstep, 1)
        m <= Setting(:forecast_block_size, 5)
        usual_model_forecast(m, :full, :none, [:histobs, :histpseudo, :histstates, :forecastobs, :forecastpseudo, :forecaststates, :forecast4qobs, :forecast4qpseudo], mb_matrix = true, check_empty_columns = false)
    end

    @testset "Ensure writing forecast centric packet runs without deprecation" begin
        write_forecast_centric_model_packet(m, :mode, :none, sections = [:estimation, :forecast, :irf])
        write_standard_model_packet(m, :mode, :none, sections = [:estimation, :forecast, :irf])
        @test_broken plot_standard_model_packet(m, :mode, :none, sections = [:estimation, :forecast, :irf])
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :forecastobs)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :bddforecastobs)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :forecastpseudo)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :forecaststates)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :shockdecobs)
        m <= Setting(:date_forecast_end, DSGE.quartertodate("2020-Q1"))
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :forecastobs)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :forecastpseudo)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :bddforecastpseudo)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :forecaststates)
        @test_broken DSGE.make_forecast_plots(m, :mode, :none, :shockdecobs)
        @test_throws ErrorException DSGE.make_forecast_plots(m, :mode, :none, :y_t)

        DSGE.plot_irf_section(m, :mode, :none, [:hist_obs])

        @test DSGE.print_variable_means(m, :none, :histobs, :obs_gdp, ["a", "b"], [quartertodate("2007-Q1")], true) == "a                                 & 0.9 \\\\\nb                                 &     \\\\\n\\end{tabular}"

        DSGE.packet_help()

        m <= Setting(:date_forecast_start, quartertodate("2019-Q1"))
        @test DSGE.month_label(m) == "Oct"
    end
end
