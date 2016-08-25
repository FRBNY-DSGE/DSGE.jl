using DSGE, Base.Test, HDF5
include("../util.jl")

# Initialize model object
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :date_conditional_end => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :forecast_kill_shocks => Setting(:forecast_kill_shocks, true))
m = Model990(custom_settings = custom_settings)
m.testing = true
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

    # 1. Test unconditional forecast

    ## Historical states and shocks
    histstates = forecast_outputs[(:none, input_type)][:histstates][1]
    histshocks = forecast_outputs[(:none, input_type)][:histshocks][1]

    df = forecast_outputs[(:none, input_type)][:df]
    data = df_to_matrix(m, df; cond_type = :none)
    sys  = compute_system(m)
    exp_histstates, exp_histshocks, _, kal = DSGE.filterandsmooth(m, data, sys)

    @test_matrix_approx_eq exp_histstates histstates
    @test_matrix_approx_eq exp_histshocks histshocks

    ## Forecasted states, observables, pseudos, and shocks
    forecaststates = forecast_outputs[(:none, input_type)][:forecaststates][1]
    forecastobs    = forecast_outputs[(:none, input_type)][:forecastobs][1]
    forecastpseudo = forecast_outputs[(:none, input_type)][:forecastpseudo][1]
    forecastshocks = forecast_outputs[(:none, input_type)][:forecastshocks][1]

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

    @test_matrix_approx_eq exp_forecaststates forecaststates
    @test_matrix_approx_eq exp_forecastobs    forecastobs
    @test_matrix_approx_eq exp_forecastpseudo forecastpseudo
    @test_matrix_approx_eq exp_forecastshocks forecastshocks

    # 2. Test conditional and semiconditional forecasts
    for cond_type in [:semi, :full]

        ## Historical states and shocks
        histstates = forecast_outputs[(cond_type, input_type)][:histstates][1]
        histshocks = forecast_outputs[(cond_type, input_type)][:histshocks][1]

        df = forecast_outputs[(cond_type, input_type)][:df]
        data = df_to_matrix(m, df; cond_type = cond_type)
        sys  = compute_system(m)
        exp_histstates, exp_histshocks, exp_histpseudo, kal = DSGE.filterandsmooth(m, data, sys)
        T = DSGE.subtract_quarters(date_forecast_start(m), date_prezlb_start(m))

        @test_matrix_approx_eq exp_histstates[:, 1:T] histstates
        @test_matrix_approx_eq exp_histshocks[:, 1:T] histshocks

        ## Forecasted states, observables, pseudos, and shocks
        forecaststates = forecast_outputs[(cond_type, input_type)][:forecaststates][1]
        forecastobs    = forecast_outputs[(cond_type, input_type)][:forecastobs][1]
        forecastpseudo = forecast_outputs[(cond_type, input_type)][:forecastpseudo][1]
        forecastshocks = forecast_outputs[(cond_type, input_type)][:forecastshocks][1]

        zend = kal[:zend]
        forecast = compute_forecast(sys[:TTT], sys[:RRR], sys[:CCC], sys[:ZZ], sys[:DD], Z_pseudo, D_pseudo,
                                    forecast_horizons(m), shocks, zend)
        exp_forecaststates = hcat(exp_histstates[:, T+1:end], forecast[:states])
        exp_forecastobs    = hcat(data[:, T+1:end],           forecast[:observables])
        exp_forecastpseudo = hcat(exp_histpseudo[:, T+1:end], forecast[:pseudo_observables])
        exp_forecastshocks = hcat(exp_histshocks[:, T+1:end], forecast[:shocks])

        @test_matrix_approx_eq exp_forecaststates forecaststates
        @test_matrix_approx_eq exp_forecastobs    forecastobs
        @test_matrix_approx_eq exp_forecastpseudo forecastpseudo
        @test_matrix_approx_eq exp_forecastshocks forecastshocks

    end # cond_type

end # input_type

nothing
