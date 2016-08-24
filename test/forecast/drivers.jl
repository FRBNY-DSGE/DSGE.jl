using DSGE, Base.Test, HDF5
include("../util.jl")

# Initialize model object
custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :date_conditional_end => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :forecast_kill_shocks => Setting(:forecast_kill_shocks, true))
m = Model990(custom_settings = custom_settings)
m.testing = true

# Run forecasts
forecast_outputs = Dict{Tuple{Symbol, Symbol}, Dict{Symbol, Any}}()
for cond_type in [:none, :semi, :full]
    df = load_data(m; cond_type=cond_type, try_disk=true, verbose=:none)
    for output_type in [:states, :shocks, :forecast]
        forecast_outputs[(cond_type, output_type)] = 
            forecast_one(m, df; input_type = :init, output_type = output_type, cond_type = cond_type)
        (forecast_outputs[(cond_type, output_type)])[:df] = df
    end
end


#############################
# Test unconditional forecast
#############################

## Historical states and shocks
histstates = forecast_outputs[(:none, :states)][:histstates][1]
histshocks = forecast_outputs[(:none, :shocks)][:histshocks][1]

df = forecast_outputs[(:none, :states)][:df]
data = df_to_matrix(m, df; cond_type = :none)
sys  = compute_system(m)
exp_histstates, exp_histshocks, _ = DSGE.filterandsmooth(m, data, sys)

@test_matrix_approx_eq exp_histstates histstates
@test_matrix_approx_eq exp_histshocks histshocks

## Forecasted states, observables, pseudos, and shocks
forecaststates = forecast_outputs[(:none, :forecast)][:forecaststates][1]
forecastobs    = forecast_outputs[(:none, :forecast)][:forecastobs][1]
forecastpseudo = forecast_outputs[(:none, :forecast)][:forecastpseudo][1]
forecastshocks = forecast_outputs[(:none, :forecast)][:forecastshocks][1]

kalman = DSGE.filter(m, data, sys)
zend = kalman[:zend]
shocks = zeros(Float64, n_shocks_exogenous(m), forecast_horizons(m))
Z_pseudo = zeros(Float64, 12, n_states_augmented(m))
D_pseudo = zeros(Float64, 12)
forecast = compute_forecast(sys[:TTT], sys[:RRR], sys[:CCC], sys[:ZZ], sys[:DD], Z_pseudo, D_pseudo, forecast_horizons(m), shocks, zend)
exp_forecaststates = forecast[:states]
exp_forecastobs    = forecast[:observables]
exp_forecastpseudo = forecast[:pseudo_observables]
exp_forecastshocks = forecast[:shocks]

@test_matrix_approx_eq exp_forecaststates forecaststates
@test_matrix_approx_eq exp_forecastobs    forecastobs
@test_matrix_approx_eq exp_forecastpseudo forecastpseudo
@test_matrix_approx_eq exp_forecastshocks forecastshocks


nothing
