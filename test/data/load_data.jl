using DSGE
using Test, DataFrames, HDF5

path = dirname(@__FILE__)

@testset "Test various forms of data loading" begin
    # Can we actually test? Require that FRED API key exists
    if haskey(ENV, "FRED_API_KEY") || isfile(joinpath(ENV["HOME"],".freddatarc"))

        # Specify vintage and dates
        global custom_settings = Dict{Symbol, Setting}(
                :data_vintage            => Setting(:data_vintage, "160812"),
                :cond_vintage            => Setting(:cond_vintage, "160812"),
                :cond_id                 => Setting(:cond_id, 0),
                :use_population_forecast => Setting(:use_population_forecast, true),
                :date_forecast_start     => Setting(:date_forecast_start, quartertodate("2016-Q3")),
                :date_conditional_end    => Setting(:date_forecast_start, quartertodate("2016-Q3")),
                :n_anticipated_shocks    => Setting(:n_anticipated_shocks, 6))

        global m = Model990(custom_settings = custom_settings, testing = true)

        # Read expected results
        exp_data, exp_cond_data, exp_semicond_data =
            jldopen("$path/../reference/load_data_out.jld2", "r") do file
                read(file, "data"), read(file, "cond_data"), read(file, "semi_cond_data")
            end

        # Unconditional data
        println("The following warnings are expected test behavior:")
        df = load_data(m; try_disk=false, verbose=:none)
        data = df_to_matrix(m, df)
        @test @test_matrix_approx_eq exp_data data

        # Conditional data
        cond_df = load_data(m; cond_type=:full, try_disk=false, verbose=:none)
        cond_data = df_to_matrix(m, cond_df; cond_type=:full)
        @test @test_matrix_approx_eq exp_cond_data cond_data

        # Semiconditional data
        semicond_df = load_data(m; cond_type=:semi, try_disk=false, verbose=:none)
        semicond_data = df_to_matrix(m, semicond_df; cond_type=:semi)
        @test @test_matrix_approx_eq exp_semicond_data semicond_data
    else
        @warn "Skipping load_data test because FRED_API_KEY not present"
    end
end

nothing
