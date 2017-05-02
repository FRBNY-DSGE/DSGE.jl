using DSGE, Base.Test, JLD

# Initialize model object
m = AnSchorfheide(testing = true)
dir = joinpath(saveroot(m), "output_data", "an_schorfheide", "ss0")

# get_forecast_input_file
@test get_forecast_input_file(m, :init) == ""
@test get_forecast_input_file(m, :mode) == rawpath(m, "estimate", "paramsmode.h5")
@test get_forecast_input_file(m, :full) == rawpath(m, "estimate", "mhsave.h5")
overrides = forecast_input_file_overrides(m)
overrides[:mode] = tempname()
@test_throws ErrorException get_forecast_input_file(m, :mode)
overrides[:full] = mktemp()[1]
@test get_forecast_input_file(m, :full) == overrides[:full]
@test get_forecast_input_file(m, :subset) == overrides[:full]
overrides[:subset] = mktemp()[1]
@test get_forecast_input_file(m, :subset) == overrides[:subset]

# get_forecast_filename
@test get_forecast_filename(m, :mode, :none, :histstates) == joinpath(dir, "forecast", "raw", "histstates_cond=none_para=mode_test.jld")
@test get_forecast_filename(m, :mode, :none, :histstates; pathfcn = workpath) == joinpath(dir, "forecast", "work", "histstates_cond=none_para=mode_test.jld")
@test get_forecast_filename(m, :mode, :none, :histstates; forecast_string = "test") == joinpath(dir, "forecast", "raw", "histstates_cond=none_fcid=test_para=mode_test.jld")
@test get_forecast_filename(m, :mode, :none, :histstates; fileformat = "h5") == joinpath(dir, "forecast", "raw", "histstates_cond=none_para=mode_test.h5")

# get_forecast_filestring_addl
@test DSGE.get_forecast_filestring_addl(:mode, :none) == ["para=mode", "cond=none"]
@test DSGE.get_forecast_filestring_addl(:mode, :none; forecast_string = "test") == ["para=mode", "cond=none", "fcid=test"]

# get_forecast_output_files
dict = Dict{Symbol, String}()
output_vars = [:histstates, :histobs, :histpseudo, :histshocks,
               :forecastobs, :shockdecobs, :dettrendobs, :trendobs, :irfobs]
for var in output_vars
    dict[var] = get_forecast_filename(m, :mode, :none, var)
end
@test get_forecast_output_files(m, :mode, :none, output_vars) == dict

# write_forecast_metadata
for var in output_vars
    jldopen(dict[var], "w") do file
        DSGE.write_forecast_metadata(m, file, var)
    end
end
jldopen(dict[:histstates], "r") do file
    dates = read(file, "date_indices")
    @test dates[date_mainsample_start(m)] == 1
    @test dates[date_mainsample_end(m)] == length(dates)
    @test read(file, "state_indices") == merge(m.endogenous_states, m.endogenous_states_augmented)
    @test all(x -> x == Symbol("DSGE.identity"), values(read(file, "state_revtransforms")))
end
jldopen(dict[:histobs], "r") do file
    @test read(file, "observable_indices") == m.observables
    @test !all(x -> x == Symbol("DSGE.identity"), values(read(file, "observable_revtransforms")))
end
jldopen(dict[:histpseudo], "r") do file
    @test read(file, "pseudoobservable_indices") == pseudo_measurement(m)[2].inds
    @test !all(x -> x == Symbol("DSGE.identity"), values(read(file, "pseudoobservable_revtransforms")))
end
jldopen(dict[:histshocks], "r") do file
    @test read(file, "shock_indices") == m.exogenous_shocks
    @test all(x -> x == Symbol("DSGE.identity"), values(read(file, "shock_revtransforms")))
end
jldopen(dict[:histshocks], "r") do file
    @test read(file, "shock_indices") == m.exogenous_shocks
    @test all(x -> x == Symbol("DSGE.identity"), values(read(file, "shock_revtransforms")))
end
jldopen(dict[:forecastobs], "r") do file
    dates = read(file, "date_indices")
    @test dates[date_forecast_start(m)] == 1
    @test dates[date_forecast_end(m)] == length(dates)
end
for class in [:shockdec, :dettrend, :trend]
    jldopen(dict[Symbol(class, "obs")], "r") do file
        dates = read(file, "date_indices")
        @test dates[date_mainsample_start(m)] == 1
        @test dates[date_forecast_end(m)] == length(dates)
    end
end
jldopen(dict[:irfobs], "r") do file
    @test !HDF5.exists(file, "date_indices")
    @test all(x -> x == Symbol("DSGE.identity"), values(read(file, "observable_revtransforms")))
end