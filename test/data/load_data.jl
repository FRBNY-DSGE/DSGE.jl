using DSGE
using Test, DataFrames, HDF5, JLD2

path = dirname(@__FILE__)

@testset "Test various forms of data loading" begin
    # Can we actually test? Require that FRED API key exists
    homedirpath = Sys.iswindows() ? joinpath(homedir(),".freddatarc") : joinpath(ENV["HOME"],".freddatarc")
    if haskey(ENV, "FRED_API_KEY") || isfile(homedirpath)

        # Specify vintage and dates
        global custom_settings = Dict{Symbol, Setting}(
                :data_vintage             => Setting(:data_vintage, "160812"),
                :cond_vintage             => Setting(:cond_vintage, "160812"),
                :cond_id                  => Setting(:cond_id, 0),
                :use_population_forecast  => Setting(:use_population_forecast, true),
                :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
                :date_conditional_end     => Setting(:date_conditional_end, DSGE.quartertodate("2016-Q3")),
                :n_mon_anticipated_shocks => Setting(:n_mon_anticipated_shocks, 6))

        global m = Model990(custom_settings = custom_settings, testing = true)
        m <= Setting(:rate_expectations_source, :ois)

        # Read expected results
        exp_data, exp_cond_data, exp_semicond_data =
            JLD2.jldopen("$path/../reference/load_data_out.jld2", "r") do file
                read(file, "data"), read(file, "cond_data"), read(file, "semi_cond_data")
            end

        # Check high summary statistics runs without an error
        @info "The following summary statistics are expected."
        load_data(m; try_disk = false, summary_statistics = :high, verbose=:none,
                  check_empty_columns = false)
        load_data(m; try_disk = false, summary_statistics = :low, verbose=:none,
                  check_empty_columns = false)

        # Unconditional data
        println("The following warnings are expected test behavior:")
        global df = load_data(m; try_disk=false, verbose=:none, summary_statistics=:none,
                              check_empty_columns = false)
        global data = df_to_matrix(m, df)
        @test @test_matrix_approx_eq exp_data data

        # Conditional data
        cond_df = load_data(m; cond_type=:full, try_disk=false, verbose=:none,
                            summary_statistics=:none,
                            check_empty_columns = false)
        cond_data = df_to_matrix(m, cond_df; cond_type=:full)
        @test @test_matrix_approx_eq exp_cond_data cond_data

        # Semiconditional data
        semicond_df = load_data(m; cond_type=:semi, try_disk=false, verbose=:none,
                                summary_statistics=:none, check_empty_columns = false)
        semicond_data = df_to_matrix(m, semicond_df; cond_type=:semi)
        @test @test_matrix_approx_eq exp_semicond_data semicond_data

        # Test errors
        m <= Setting(:data_vintage, "160813")
        m <= Setting(:cond_vintage, "160812")
        @test_throws ErrorException load_cond_data_levels(m)
        m <= Setting(:cond_vintage, "160813")
        @test_throws ErrorException load_cond_data_levels(m)
        @test_throws ErrorException load_data(m; try_disk = false, verbose=:none)
    else
        @warn "Skipping load_data test because FRED_API_KEY not present"
    end
end

nothing
