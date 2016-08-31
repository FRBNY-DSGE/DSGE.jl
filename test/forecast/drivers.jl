using DSGE, Base.Test, HDF5
include("../util.jl")

# Initialize model object
custom_settings = Dict{Symbol, Setting}(
    :data_vintage            => Setting(:data_vintage, "160812"),
    :cond_vintage            => Setting(:cond_vintage, "160812"),
    :use_population_forecast => Setting(:use_population_forecast, true),
    :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2016-Q3")),
    :date_conditional_end    => Setting(:date_conditional_end, quartertodate("2016-Q3")),
    :date_forecast_end       => Setting(:date_forecast_end, quartertodate("2016-Q4")),
    :n_anticipated_shocks    => Setting(:n_anticipated_shocks, 6),
    # :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    # :date_conditional_end => Setting(:date_conditional_end, quartertodate("2015-Q4")),
    :forecast_kill_shocks => Setting(:forecast_kill_shocks, true))
m = Model990(custom_settings = custom_settings, testing = true)
init_params = map(θ -> θ.value, m.parameters)

# Copy mode
mode_infile = "paramsmode_m990_ss2_nant=0_vint=151127.h5"
mode_inpath = joinpath(m.settings[:saveroot].value, "input_data/user", mode_infile)
mode_params = h5open(mode_inpath, "r") do h5
    read(h5, "params")
end

mode_outpath = DSGE.get_input_file(m, :mode)
h5open(mode_outpath, "w") do h5
    write(h5, "params", mode_params)
end

# Run forecasts
forecast_outputs = Dict{Tuple{Symbol, Symbol}, Dict{Symbol, Any}}()
output_vars = [:histstates, :histpseudo, :histshocks, :forecaststates, :forecastpseudo, :forecastobs, :forecastshocks]

for input_type in [:init, :mode]
    for cond_type in [:none, :semi, :full]

        forecast_output = Dict{Symbol, Any}()
        forecast_output[:df] = load_data(m; cond_type=cond_type, try_disk=true, verbose=:none)

        new_forecast = forecast_one(m, forecast_output[:df]; input_type =
            input_type, cond_type = cond_type, output_vars = output_vars)
        merge!(forecast_output, new_forecast)

        forecast_outputs[(cond_type, input_type)] = forecast_output
    end
end

# Test forecast outputs
for input_type in [:init, :mode]

    if input_type == :init
        update!(m, init_params)
    elseif input_type == :mode
        update!(m, mode_params)
    end

    for cond_type in [:none, :semi, :full]

        ## Historical states and shocks
        histstates = forecast_outputs[(cond_type, input_type)][:histstates][1]
        histshocks = forecast_outputs[(cond_type, input_type)][:histshocks][1]

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
        forecaststates = forecast_outputs[(cond_type, input_type)][:forecaststates][1]
        forecastobs    = forecast_outputs[(cond_type, input_type)][:forecastobs][1]
        forecastpseudo = forecast_outputs[(cond_type, input_type)][:forecastpseudo][1]
        forecastshocks = forecast_outputs[(cond_type, input_type)][:forecastshocks][1]

        zend = kal[:zend]
        shocks = zeros(Float64, n_shocks_exogenous(m), forecast_horizons(m))
        Z_pseudo = zeros(Float64, 12, n_states_augmented(m))
        D_pseudo = zeros(Float64, 12)
        forecast = compute_forecast(sys[:TTT], sys[:RRR], sys[:CCC], sys[:ZZ], sys[:DD], Z_pseudo, D_pseudo,
                                    forecast_horizons(m), shocks, zend)

        exp_forecaststates = forecast[:states]
        exp_forecastobs    = forecast[:observables]
        exp_forecastpseudo = forecast[:pseudo_observables]
        exp_forecastshocks = forecast[:shocks]

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

    end # cond_type

end # input_type

nothing
