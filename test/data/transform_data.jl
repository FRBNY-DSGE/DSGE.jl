# Load data to use for tests
path = dirname(@__FILE__)
fred = CSV.read("$path/../reference/fred_160812.csv", DataFrame)

# Specify vintage and dates
custom_settings = Dict{Symbol, Setting}(
    :data_vintage             => Setting(:data_vintage, "160812"),
    :cond_vintage             => Setting(:cond_vintage, "160812"),
    :cond_id                  => Setting(:cond_id, 0),
    :use_population_forecast  => Setting(:use_population_forecast, true),
    :date_forecast_start      => Setting(:date_forecast_start, DSGE.quartertodate("2016-Q3")),
    :date_conditional_end     => Setting(:date_conditional_end, DSGE.quartertodate("2016-Q3")),
    :n_mon_anticipated_shocks => Setting(:n_mon_anticipated_shocks, 6))
m = Model990(testing = true, custom_settings = custom_settings)
m <= Setting(:rate_expectations_source, :ois)

# Test full, semi, and none conditional data
exp_data, exp_cond_data, exp_semicond_data =
    JLD2.jldopen("$path/../reference/load_data_out.jld2", "r") do file
        read(file, "data"), read(file, "cond_data"), read(file, "semi_cond_data")
    end
levels, semi_levels, full_levels =
    JLD2.jldopen("$path/../reference/transform_data_inputs.jld2", "r") do file
        read(file, "none"), read(file, "semi"), read(file, "full")
end
exp_data_rev_transforms =
    JLD2.jldopen("$path/../reference/transform_data_out.jld2", "r") do file
        read(file, "rev_transform")
    end
exp_hp_pop_hist, exp_hp_pop_forecast =
    JLD2.jldopen("$path/../reference/hp_data_out.jld2", "r") do file
        read(file, "exp_hp_pop_hist"), read(file, "exp_hp_pop_forecast")
    end
pop_forecast = DataFrame(CNP16OV = cumsum(fred[!,:CNP16OV][end] .* (ones(8) .* 1.03)))
pop_forecast[!,:date] = DSGE.get_quarter_ends(DSGE.next_quarter(fred[!,:date][end]), fred[!,:date][end] + Year(2))

@testset "Check helper functions for transform_data" begin
    @test !isempty(collect_data_transforms(m))
    @test collect_data_transforms(m, direction = :rev) == exp_data_rev_transforms
    hp_hist, hp_forecast = DSGE.transform_population_data(fred, pop_forecast, :CNP16OV,
                                                     pad_forecast_start = true)
    @test hp_hist[!,:filtered_population_recorded] == exp_hp_pop_hist
    @test hp_forecast[!,:filtered_population_forecast] == exp_hp_pop_forecast
    hp_hist, hp_forecast = DSGE.transform_population_data(fred, pop_forecast, :CNP16OV)
    @test hp_hist[!,:filtered_population_recorded] == exp_hp_pop_hist
    @test hp_forecast[!,:filtered_population_forecast] == exp_hp_pop_forecast

    hist, forecast = DSGE.transform_population_data(fred, pop_forecast, :CNP16OV,
                                    pad_forecast_start = true, use_hpfilter = false)
    @test @test_matrix_approx_eq hist[!,:dlpopulation_recorded] difflog(Float64.(fred[!,:CNP16OV]))
    @test @test_matrix_approx_eq forecast[!,:dlpopulation_forecast] difflog(Float64.(pop_forecast[!,:CNP16OV]))
    hist, forecast = DSGE.transform_population_data(fred, pop_forecast, :CNP16OV,
                                    use_hpfilter = false)
    @test @test_matrix_approx_eq hist[!,:dlpopulation_recorded] difflog(Float64.(fred[!,:CNP16OV]))
    @test @test_matrix_approx_eq forecast[!,:dlpopulation_forecast] difflog(Float64.(pop_forecast[!,:CNP16OV]))
end

@testset "Test transform_data for different conditional data" begin
    @info "The following warnings are expected."
    df = transform_data(m, levels; cond_type = :none, verbose = :none)
    start_date = date_presample_start(m)
    end_date = date_mainsample_end(m)
    df = df[start_date .<= df[!,:date] .<= end_date, :]
    DSGE.missing_cond_vars!(m, df, cond_type = :none, check_empty_columns = false)
    @test @test_matrix_approx_eq df_to_matrix(m,df) exp_data
    df = transform_data(m, semi_levels; cond_type = :semi, verbose = :none)
    end_date = date_conditional_end(m)
    df = df[start_date .<= df[!,:date] .<= end_date, :]
    DSGE.missing_cond_vars!(m, df, cond_type = :semi, check_empty_columns = false)
    @test @test_matrix_approx_eq df_to_matrix(m,df, cond_type = :semi) exp_semicond_data
    df = transform_data(m, full_levels; cond_type = :full, verbose = :none)
    end_date = date_conditional_end(m)
    df = df[start_date .<= df[!,:date] .<= end_date, :]
    DSGE.missing_cond_vars!(m, df, cond_type = :full, check_empty_columns = false)
    @test @test_matrix_approx_eq df_to_matrix(m,df, cond_type = :full) exp_cond_data
end

nothing
