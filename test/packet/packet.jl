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
        # These plotting calls were @test_broken in 2020 because they errored then;
        # under Julia 1.12 the individual forecast/pseudo/state plots now run and
        # return an OrderedDict of plots, so those are smoke tests that they execute
        # without error. Calls that depend on output this run doesn't produce stay
        # @test_broken: plot_standard_model_packet (its :forecast section plots shock
        # decompositions, whose mbshockdec*.jld2 is never generated here) and the
        # bare :shockdecobs plots. Guarded with `; false` so a non-Boolean return, if
        # the behavior ever changes, can't crash the testset.
        @test_broken (plot_standard_model_packet(m, :mode, :none, sections = [:estimation, :forecast, :irf]); false)
        @test (DSGE.make_forecast_plots(m, :mode, :none, :forecastobs); true)
        @test (DSGE.make_forecast_plots(m, :mode, :none, :bddforecastobs); true)
        @test (DSGE.make_forecast_plots(m, :mode, :none, :forecastpseudo); true)
        @test (DSGE.make_forecast_plots(m, :mode, :none, :forecaststates); true)
        @test_broken (DSGE.make_forecast_plots(m, :mode, :none, :shockdecobs); false)
        m <= Setting(:date_forecast_end, DSGE.quartertodate("2020-Q1"))
        @test (DSGE.make_forecast_plots(m, :mode, :none, :forecastobs); true)
        @test (DSGE.make_forecast_plots(m, :mode, :none, :forecastpseudo); true)
        @test_broken (DSGE.make_forecast_plots(m, :mode, :none, :bddforecastpseudo); false)
        @test (DSGE.make_forecast_plots(m, :mode, :none, :forecaststates); true)
        @test_broken (DSGE.make_forecast_plots(m, :mode, :none, :shockdecobs); false)
        @test_throws ErrorException DSGE.make_forecast_plots(m, :mode, :none, :y_t)

        DSGE.plot_irf_section(m, :mode, :none, [:hist_obs])

        @test DSGE.print_variable_means(m, :none, :histobs, :obs_gdp, ["a", "b"], [quartertodate("2007-Q1")], true) == "a                                 & 0.9 \\\\\nb                                 &     \\\\\n\\end{tabular}"

        DSGE.packet_help()

        m <= Setting(:date_forecast_start, quartertodate("2019-Q1"))
        @test DSGE.month_label(m) == "Oct"
    end
end
