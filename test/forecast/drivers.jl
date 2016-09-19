using DSGE, Base.Test, HDF5
include("../util.jl")

# Initialize model object
custom_settings = Dict{Symbol, Setting}(
    :use_population_forecast => Setting(:use_population_forecast, true),
    :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :date_conditional_end    => Setting(:date_conditional_end, quartertodate("2015-Q4")),
    :date_forecast_end       => Setting(:date_forecast_end, quartertodate("2016-Q1")),
    :n_anticipated_shocks    => Setting(:n_anticipated_shocks, 6),
    :forecast_kill_shocks    => Setting(:forecast_kill_shocks, true),
    :saveroot                => Setting(:saveroot, normpath(joinpath(dirname(@__FILE__), "..", "reference"))),
    :use_parallel_workers    => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings, testing = true)

# Add parallel workers
my_procs = addprocs(10)
@everywhere using DSGE

# Run forecasts
forecast_outputs = Dict{Tuple{Symbol, Symbol}, Dict{Symbol, Any}}()
output_vars = [:histstates, :histpseudo, :histshocks, :forecaststates,
               :forecastpseudo, :forecastobs, :forecastshocks, :shockdecstates,
               :shockdecpseudo, :shockdecobs]
output_files = []

for input_type in [:init, :mode, :full]

    # Call forecast_one once without timing
    df = load_data(m; verbose = :none)
    forecast_one(m, df; input_type = :full, cond_type = :none, output_vars = output_vars)

    for cond_type in [:none, :semi, :full]

        println("input_type = $(input_type), cond_type = $(cond_type)")

        forecast_output = Dict{Symbol, Any}()
        forecast_output[:df] = load_data(m; cond_type=cond_type, try_disk=true, verbose=:none)

        @time new_forecast = forecast_one(m, forecast_output[:df]; input_type =
            input_type, cond_type = cond_type, output_vars = output_vars)
        merge!(forecast_output, new_forecast)

        forecast_outputs[(cond_type, input_type)] = forecast_output

        new_output_files = collect(values(DSGE.get_output_files(m, input_type, output_vars, cond_type)))
        append!(output_files, new_output_files)
    end
end

# Test forecast outputs
for input_type in [:init, :mode]

    if input_type == :init
        DSGE.init_parameters!(m)
        DSGE.steadystate!(m)
    elseif input_type == :mode
        specify_mode!(m, DSGE.get_input_file(m, :mode))
    end

    for cond_type in [:none, :semi, :full]

        ## Historical states and shocks
        histstates = forecast_outputs[(cond_type, input_type)][:histstates]
        histshocks = forecast_outputs[(cond_type, input_type)][:histshocks]

        df = forecast_outputs[(cond_type, input_type)][:df]
        data = df_to_matrix(m, df; cond_type = cond_type)
        sys  = compute_system(m)
        exp_histstates, exp_histshocks, exp_histpseudo, kal = DSGE.filterandsmooth(m, data, sys)

        if cond_type in [:semi, :full]
            T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))
            @test_matrix_approx_eq exp_histstates[:, 1:T] histstates
            @test_matrix_approx_eq exp_histshocks[:, 1:T] histshocks
        else
            @test_matrix_approx_eq exp_histstates histstates
            @test_matrix_approx_eq exp_histshocks histshocks
        end            

        ## Forecasted states, observables, pseudos, and shocks
        forecaststates = forecast_outputs[(cond_type, input_type)][:forecaststates]
        forecastobs    = forecast_outputs[(cond_type, input_type)][:forecastobs]
        forecastpseudo = forecast_outputs[(cond_type, input_type)][:forecastpseudo]
        forecastshocks = forecast_outputs[(cond_type, input_type)][:forecastshocks]

        zend = kal[:zend]
		shocks = zeros(Float64, n_shocks_exogenous(m), forecast_horizons(m))
        _, pseudo_measur = pseudo_measurement(m)
        Z_pseudo = pseudo_measur.ZZ
        D_pseudo = pseudo_measur.DD        
        exp_forecaststates, exp_forecastobs, exp_forecastpseudo, exp_forecastshocks =
            DSGE.compute_forecast(sys[:TTT], sys[:RRR], sys[:CCC], sys[:ZZ], sys[:DD],
                                  Z_pseudo, D_pseudo, forecast_horizons(m), shocks, zend)

        if cond_type in [:semi, :full]
            exp_histobs = data[:, index_prezlb_start(m):end]
            @test_matrix_approx_eq hcat(exp_histstates[:, T+1:end], exp_forecaststates) forecaststates
            @test_matrix_approx_eq hcat(exp_histobs[:, T+1:end],    exp_forecastobs)    forecastobs
            @test_matrix_approx_eq hcat(exp_histpseudo[:, T+1:end], exp_forecastpseudo) forecastpseudo
            @test_matrix_approx_eq hcat(exp_histshocks[:, T+1:end], exp_forecastshocks) forecastshocks
        else
            @test_matrix_approx_eq exp_forecaststates forecaststates
            @test_matrix_approx_eq exp_forecastobs    forecastobs
            @test_matrix_approx_eq exp_forecastpseudo forecastpseudo
            @test_matrix_approx_eq exp_forecastshocks forecastshocks
        end

        ## Shock decompositions of states, pseudo-observables, and observables
        shockdecstates = forecast_outputs[(cond_type, input_type)][:shockdecstates]
        shockdecobs    = forecast_outputs[(cond_type, input_type)][:shockdecobs]
        shockdecpseudo = forecast_outputs[(cond_type, input_type)][:shockdecpseudo]

        exp_shockdecstates, exp_shockdecobs, exp_shockdecpseudo =
            DSGE.compute_shock_decompositions(sys[:TTT], sys[:RRR], sys[:ZZ],
                                              sys[:DD], Z_pseudo, D_pseudo, forecast_horizons(m), exp_histshocks)

        @test_matrix_approx_eq exp_shockdecstates shockdecstates
        @test_matrix_approx_eq exp_shockdecobs    shockdecobs
        @test_matrix_approx_eq exp_shockdecpseudo shockdecpseudo

    end # cond_type

end # input_type

# Delete all files written by forecast_one
map(rm, unique(output_files))

# Remove parallel workers
rmprocs(my_procs)

nothing
