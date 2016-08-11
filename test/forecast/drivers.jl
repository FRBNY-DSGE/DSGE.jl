using DSGE

custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_population_forecast => Setting(:use_population_forecast, true))
m = Model990(custom_settings = custom_settings)
m.testing = true

df = load_data(m; try_disk=true, verbose=:none)
for output_type in [:states, :shocks, :simple, :forecast]
    for cond_type in [:none, :semi, :full]
        forecast_output = forecast_one(m, df; input_type = :init, output_type = output_type, cond_type = cond_type)
    end
end

nothing
