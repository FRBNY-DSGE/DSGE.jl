using DSGE
using Base.Test, DataFrames

# Specify exact vintage and dates
custom_settings = Dict{Symbol, Setting}(
    :data_vintage            => Setting(:data_vintage, "151127"),
    :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2015-Q4")),
    :use_population_forecast => Setting(:use_population_forecast, true))
m = Model990(custom_settings = custom_settings)
m.testing = true

# Conditional data
cond_df = load_cond_data(m, :full, verbose=:high)

@assert cond_df[end, :date] == date_forecast_start(m)
cond_obs = [:obs_gdp, :obs_corepce, :obs_nominalrate, :obs_spread]
for obs in cond_obs
    @assert !isnan(cond_df[end, obs])
end
for obs in setdiff(names(cond_df), [cond_obs; :date])
    @assert isnan(cond_df[end, obs])
end

# Semiconditional data
semicond_df = load_cond_data(m, :semi; verbose=:high)

@assert semicond_df[end, :date] == date_forecast_start(m)
semicond_obs = get_setting(m, :cond_semi_names)
for obs in semicond_obs
    @assert !isnan(semicond_df[end, obs])
end
for obs in setdiff(names(semicond_df), [semicond_obs; :date])
    @assert isnan(semicond_df[end, obs])
end

nothing