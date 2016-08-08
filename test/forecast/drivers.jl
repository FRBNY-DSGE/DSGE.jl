using DSGE

custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = Model990(custom_settings = custom_settings)
m.testing = true

df = load_data(m; try_disk=true, verbose=:none)
for output_type in [:states, :shocks, :simple, :forecast]
    forecast_output = forecast_one(m, df; input_type = :init, output_type = output_type, cond_type = :none)
end

nothing
