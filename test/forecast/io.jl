using DSGE, Test, Dates
using JLD2, OrderedCollections, HDF5, Nullables

path = dirname(@__FILE__)
save_output = true

# Initialize model object
m = AnSchorfheide(testing = true)
dir = joinpath(saveroot(m), "output_data", "an_schorfheide", "ss0")

output_vars = [:histstates, :histobs, :histpseudo, :histshocks,
                   :forecastobs, :shockdecobs, :dettrendobs, :trendobs, :irfobs]

@testset "Test main forecast I/O functions and metadata writing" begin
    # get_forecast_input_file
    @test get_forecast_input_file(m, :init) == ""
    @test get_forecast_input_file(m, :init_draw_shocks) == ""
    @test get_forecast_input_file(m, :prior) == ""
    @test get_forecast_input_file(m, :mode) == rawpath(m, "estimate", "paramsmode.h5")
    @test get_forecast_input_file(m, :mode_draw_shocks) == rawpath(m, "estimate", "paramsmode.h5")
    m <= Setting(:sampling_method, :invalid_method)
    @test_throws ErrorException get_forecast_input_file(m, :full) == rawpath(m, "estimate", "smcsave.h5")
    m <= Setting(:sampling_method, :SMC)
    @test get_forecast_input_file(m, :full) == rawpath(m, "estimate", "smcsave.h5")
    m <= Setting(:sampling_method, :MH)
    @test get_forecast_input_file(m, :full) == rawpath(m, "estimate", "mhsave.h5")
    global overrides = forecast_input_file_overrides(m)
    overrides[:full] = "invalidpath"
    @test_throws ErrorException get_forecast_input_file(m, :full)
    overrides[:mode] = tempname()
    @test_throws ErrorException get_forecast_input_file(m, :mode)
    overrides[:full] = mktemp()[1]
    @test get_forecast_input_file(m, :full) == overrides[:full]
    @test get_forecast_input_file(m, :subset) == overrides[:full]
    overrides[:subset] = mktemp()[1]
    @test get_forecast_input_file(m, :subset) == overrides[:subset]

    # get_forecast_filename
    @test get_forecast_filename(m, :mode, :none, :histstates) == joinpath(dir, "forecast", "raw", "histstates_cond=none_para=mode_test.jld2")
    @test get_forecast_filename(m, :mode, :none, :histstates; pathfcn = workpath) == joinpath(dir, "forecast", "work", "histstates_cond=none_para=mode_test.jld2")
    @test get_forecast_filename(m, :mode, :none, :histstates; forecast_string = "test") == joinpath(dir, "forecast", "raw", "histstates_cond=none_fcid=test_para=mode_test.jld2")
    @test get_forecast_filename(m, :mode, :none, :histstates; fileformat = :h5) == joinpath(dir, "forecast", "raw", "histstates_cond=none_para=mode_test.h5")

    # get_forecast_filestring_addl
    @test DSGE.get_forecast_filestring_addl(:mode, :none) == ["para=mode", "cond=none"]
    @test DSGE.get_forecast_filestring_addl(:mode, :none; forecast_string = "test") == ["para=mode", "cond=none", "fcid=test"]
    @test_throws ErrorException DSGE.get_forecast_filestring_addl(:subset, :none)

    # get_forecast_output_files
    dict = Dict{Symbol, String}()

    for var in output_vars
        dict[var] = get_forecast_filename(m, :mode, :none, var)
    end
    @test get_forecast_output_files(m, :mode, :none, output_vars) == dict

    # write_forecast_outputs
    global df = load_data(m, verbose = :none)
    overrides = forecast_input_file_overrides(m)
    overrides[:full] = "$(path)/../reference/smcsave_.h5"
    m <= Setting(:sampling_method, :SMC)
    m <= Setting(:forecast_block_size, 100)
    block_inds, block_inds_thin = DSGE.forecast_block_inds(m, :full)
    forecast_output_files = DSGE.get_forecast_output_files(m, :full, :none, [:histstates, :histobs, :forecaststates, :forecastobs])
    # Test that write_forecast_outputs and combine_raw_forecast_output_and_metadata run without deprecration
    for block in 1:20
        forecast_output = load("$(path)/../reference/write_forecast_outputs_$(block)_in.jld2", "forecast_output")
        DSGE.write_forecast_outputs(m, :full, [:histstates, :histobs, :forecaststates, :forecastobs],forecast_output_files, forecast_output, df = df, block_number = Nullables.Nullable(block), block_inds = block_inds_thin[block])
    end
    DSGE.combine_raw_forecast_output_and_metadata(m, forecast_output_files, verbose = :none)

    # write_forecast_metadata
    for var in output_vars
        JLD2.jldopen(dict[var], "w") do file
            DSGE.write_forecast_metadata(m, file, var)
        end
    end
    JLD2.jldopen(dict[:histstates], "r") do file
        dates = read(file, "date_indices")
        @test dates[date_mainsample_start(m)] == 1
        @test dates[date_mainsample_end(m)] == length(dates)
        @test read(file, "state_indices") == merge(m.endogenous_states, m.endogenous_states_augmented)
        @test all(x -> x == Symbol("identity"), values(read(file, "state_revtransforms")))
    end
    JLD2.jldopen(dict[:histobs], "r") do file
        @test read(file, "observable_indices") == m.observables
        @test !all(x -> x == Symbol("identity"), values(read(file, "observable_revtransforms")))
    end
    JLD2.jldopen(dict[:histpseudo], "r") do file
        @test read(file, "pseudoobservable_indices") == m.pseudo_observables
        @test !all(x -> x == Symbol("identity"), values(read(file, "pseudoobservable_revtransforms")))
    end
    JLD2.jldopen(dict[:histshocks], "r") do file
        @test read(file, "shock_indices") == m.exogenous_shocks
        @test all(x -> x == Symbol("identity"), values(read(file, "shock_revtransforms")))
    end
    JLD2.jldopen(dict[:histshocks], "r") do file
        @test read(file, "shock_indices") == m.exogenous_shocks
        @test all(x -> x == Symbol("identity"), values(read(file, "shock_revtransforms")))
    end
    JLD2.jldopen(dict[:forecastobs], "r") do file
        dates = read(file, "date_indices")
        @test dates[date_forecast_start(m)] == 1
        @test dates[date_forecast_end(m)] == length(dates)
    end
    for class in [:shockdec, :dettrend, :trend]
        JLD2.jldopen(dict[Symbol(class, "obs")], "r") do file
            dates = read(file, "date_indices")
            @test dates[date_mainsample_start(m)] == 1
            @test dates[date_forecast_end(m)] == length(dates)
        end
    end
    JLD2.jldopen(dict[:irfobs], "r") do file
        @test !haskey(file, "date_indices")
        @test all(x -> x == Symbol("identity"), values(read(file, "observable_revtransforms")))
    end
end
