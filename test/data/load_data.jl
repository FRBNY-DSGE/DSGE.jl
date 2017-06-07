using DSGE
using Base.Test, DataFrames, HDF5

path = dirname(@__FILE__)

# Can we actually test? Require that FRED API key exists
if haskey(ENV, "FRED_API_KEY") || isfile(joinpath(ENV["HOME"],".freddatarc"))

    # Specify vintage and dates
    custom_settings = Dict{Symbol, Setting}(
            :data_vintage            => Setting(:data_vintage, "160812"),
            :cond_vintage            => Setting(:cond_vintage, "160812"),
            :use_population_forecast => Setting(:use_population_forecast, true),
            :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2016-Q3")),
            :date_conditional_end    => Setting(:date_forecast_start, quartertodate("2016-Q3")),
            :n_anticipated_shocks    => Setting(:n_anticipated_shocks, 6))

    m = Model990(custom_settings = custom_settings, testing = true)

    # Read expected results
    exp_data, exp_cond_data, exp_semicond_data =
        h5open("$path/../reference/load_data_out.h5", "r") do h5
            read(h5, "data"), read(h5, "cond_data"), read(h5, "semicond_data")
        end

    # Unconditional data
    println("The following warnings are expected test behavior:")
    df = load_data(m; try_disk=false, verbose=:none)
    data = df_to_matrix(m, df)
    @test_matrix_approx_eq exp_data data

    # Conditional data
    cond_df = load_data(m; cond_type=:full, try_disk=false, verbose=:none)
    cond_data = df_to_matrix(m, cond_df; cond_type=:full)
    @test_matrix_approx_eq exp_cond_data cond_data

    # Semiconditional data
    semicond_df = load_data(m; cond_type=:semi, try_disk=false, verbose=:none)
    semicond_data = df_to_matrix(m, semicond_df; cond_type=:semi)
    @test_matrix_approx_eq exp_semicond_data semicond_data

    # include_presample flag
    data = df_to_matrix(m, df; include_presample = false)
    @test_matrix_approx_eq exp_data[:, 3:end] data
else
    warn("Skipping load_data test because FRED_API_KEY not present")
end
nothing