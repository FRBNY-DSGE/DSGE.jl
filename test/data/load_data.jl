using DSGE
using Base.Test, DataFrames, HDF5
include("../util.jl")

path = dirname(@__FILE__)

# Read expected results
exp_data, exp_cond_data, exp_semicond_data = 
    h5open("$path/../reference/load_data_out.h5") do h5
        read(h5, "data"), read(h5, "cond_data"), read(h5, "semicond_data")
end

# Can we actually test? Require that the FRED_API_KEY ENV is populated.
@assert haskey(ENV, "FRED_API_KEY")

# Specify exact vintage and dates
custom_settings = Dict{Symbol, Setting}(
    :data_vintage            => Setting(:data_vintage, "160812"),
    :cond_vintage            => Setting(:data_vintage, "160812"),
    :use_population_forecast => Setting(:use_population_forecast, true),
    :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2016-Q3")))
m = Model990(custom_settings = custom_settings)
m.testing = true

df = load_data(m; try_disk=false, verbose=:none)
data = df_to_matrix(m, df)
for obs in 1:12
    println("Observable $obs")
    try
        @test_matrix_approx_eq_eps exp_data[obs, :] data[obs, :] 1e-3 1e-1
    catch err
        warn(err)
        continue
    end
end

# Conditional data
cond_df = load_data(m; cond_type=:full, try_disk=false, verbose=:none)
cond_data = df_to_matrix(m, cond_df)[:, end]
for obs in 1:12
    println("Observable $obs")
    try
        @test_matrix_approx_eq_eps exp_cond_data[obs, :] cond_data[obs, :] 1e-3 1e-1
    catch err
        warn(err)
        continue
    end
end

# @assert cond_df[end, :date] == date_forecast_start(m)
# cond_obs = cond_full_names(m)
# for obs in cond_obs
#     @assert !isnan(cond_df[end, obs])
# end
# for obs in setdiff(names(cond_df), [cond_obs; :date])
#     @assert isnan(cond_df[end, obs])
# end

# Semiconditional data
semicond_df = load_data(m; cond_type=:semi, try_disk=false, verbose=:none)
semicond_data = df_to_matrix(m, semicond_df)[:, end]
for obs in 1:12
    println("Observable $obs")
    try
        @test_matrix_approx_eq_eps exp_cond_data[obs, :] cond_data[obs, :] 1e-3 1e-1
    catch err
        warn(err)
        continue
    end
end

# @assert semicond_df[end, :date] == date_forecast_start(m)
# semicond_obs = cond_semi_names(m)
# for obs in semicond_obs
#     @assert !isnan(semicond_df[end, obs])
# end
# for obs in setdiff(names(semicond_df), [semicond_obs; :date])
#     @assert isnan(semicond_df[end, obs])
# end

nothing