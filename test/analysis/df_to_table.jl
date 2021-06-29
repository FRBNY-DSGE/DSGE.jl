save_output = false

using DSGE, ModelConstructors, Test, JLD2, FileIO, OrderedCollections, CSV, Random, DataFrames
path = dirname(@__FILE__)

m = AnSchorfheide()

m <= Setting(:use_population_forecast, false)
m <= Setting(:cond_id, 0)
m <= Setting(:date_forecast_start, quartertodate("2015-Q4"))
m <= Setting(:date_conditional_end, quartertodate("2015-Q4"))
m <= Setting(:forecast_block_size, 5)
m <= Setting(:data_vintage, "REF", true, "vint", "Date of data")
m <= Setting(:dataroot, "$path/../reference/input_data")

estroot = normpath(joinpath(dirname(@__FILE__), "..", "reference"))
overrides = forecast_input_file_overrides(m)
overrides[:mode] = joinpath(estroot, "optimize.h5")
overrides[:full] = joinpath(estroot, "metropolis_hastings.h5")

output_vars = add_requisite_output_vars([:histpseudo, :histobs, :histstdshocks,
                                         :histutpseudo, :histutobs,
                                         :hist4qpseudo, :hist4qobs,
                                         :forecaststates, :forecastpseudo, :forecastobs, :forecaststdshocks,
                                         :forecastutpseudo, :forecastutobs,
                                         :forecast4qpseudo, :forecast4qobs,
                                         :bddforecaststates, :bddforecastshocks, :bddforecastpseudo, :bddforecastobs,
                                         :shockdecpseudo, :shockdecobs,
                                         :trendstates, :trendobs, :trendpseudo,
                                         :dettrendstates, :dettrendobs, :dettrendpseudo,
                                         :irfstates, :irfpseudo, :irfobs])

Random.seed!(47)
df = CSV.read("$path/../reference/input_data/data/data_dsid=00_vint=REF.csv", DataFrame)
forecast_one(m, :full, :none, output_vars, df = df, verbose = :none)
compute_meansbands(m, :full, :none, [:histobs, :forecastobs, :histpseudo, :forecastpseudo,
                                     :hist4qobs, :forecast4qobs, :hist4qpseudo, :forecast4qpseudo, :shockdecobs, :irfobs], df = df, verbose = :none)

hist, fore = construct_fcast_and_hist_dfs(m, :none, [:obs_gdp, :obs_cpi, :obs_nominalrate])
hist_noT, fore_noT = construct_fcast_and_hist_dfs(m, :none, [:obs_gdp, :obs_cpi, :obs_nominalrate], include_T_in_df_forecast = false)
hist_4q, fore_4q = construct_fcast_and_hist_dfs(m, :none, [:obs_gdp, :obs_cpi, :obs_nominalrate], use_4q = true)


if save_output
    fn = VERSION >= v"1.6" ? "$path/../reference/df_to_table_out_v1p6.jld2" : VERSION >= v"1.5" ? "$path/../reference/df_to_table_out_v1p5.jld2" : "$path/../reference/df_to_table_out.jld2"
    jldopen(fn, true, true, true, IOStream) do file
        file["hist"] = hist
        file["hist_noT"] = hist_noT
        file["hist_4q"] = hist_4q
        file["fore"] = fore
        file["fore_noT"] = fore_noT
        file["fore_4q"] = fore_4q
    end
end


fn = VERSION >= v"1.6" ? "$path/../reference/df_to_table_out_v1p6.jld2" : VERSION >= v"1.5" ? "$path/../reference/df_to_table_out_v1p5.jld2" : "$path/../reference/df_to_table_out.jld2"
saved_hist = load(fn, "hist")
saved_hist_noT = load(fn, "hist_noT")
saved_hist_4q = load(fn, "hist_4q")
saved_fore = load(fn, "fore")
saved_fore_noT = load(fn, "fore_noT")
saved_fore_4q = load(fn, "fore_4q")

@testset "Test df_to_table" begin
    @test Matrix(hist) == Matrix(saved_hist)
    @test Matrix(hist_noT) == Matrix(saved_hist_noT)
    @test Matrix(hist_4q) == Matrix(saved_hist_4q)
    @test Matrix(fore) == Matrix(saved_fore)
    @test Matrix(fore_noT) == Matrix(saved_fore_noT)
    @test Matrix(fore_4q) == Matrix(saved_fore_4q)
    @test_throws AssertionError construct_fcast_and_hist_dfs(m, :none, [:obs_gdp], save_to_table = true)
    # Test that saving to tex table functionality runs
    construct_fcast_and_hist_dfs(m, :none, [:obs_gdp], save_to_table = true, table_caption = "Test Caption",
                                 filename = "test.tex", savedir = "$path")
end

@testset "Test auxiliary methods" begin
    header_mappings = DSGE.create_table_header_mappings(m, collect(keys(m.observables)))
    @test header_mappings == Dict(:obs_gdp => Symbol("Real GDP Growth"), :obs_cpi => Symbol("CPI Inflation"), :obs_nominalrate => Symbol("Nominal FFR"))
    #df = load_data(m)
    @test DSGE.unit_mappings(m, df, header_mappings) == OrderedDict(Symbol("Real GDP Growth") =>
                                                             "(Q/Q) \\% Annualized",
                                                             Symbol("CPI Inflation") =>
                                                             "(Q/Q) \\% Annualized")
end

#These really belong in analysis/io.jl MeansBands tests but put them here to avoid duplicate meansbands computations
@testset "Test read_mb" begin
    # just make sure they run
    read_mb(m, :full, :none, :forecastobs)
    read_mb(m, :full, :none, :forecastobs, bdd_and_unbdd = true)
    DSGE.get_mb_metadata(m, :full, :none, :forecastobs)
    DSGE.get_mb_metadata(m, :full, :none, :irfobs)
end

@testset "Test creating q4q4 mb" begin
    mb_4q = read_mb(m, :full, :none, :forecast4qobs)
    @test typeof(create_q4q4_mb(mb_4q)) == MeansBands
end

@testset "Test write_meansbands_tables functions" begin
    @test_throws AssertionError write_meansbands_tables_timeseries(m, :mode, :none, :irf)
    write_meansbands_tables_timeseries(m, :full, :none, :histobs)
    write_meansbands_tables_timeseries(m, :full, :none, :histforecastobs)

    write_means_tables_shockdec(m, :full, :none, :obs)

    write_meansbands_tables_all(m, :full, :none, [:histobs, :shockdecobs])
    @test_throws ErrorException write_meansbands_tables_all(m, :full, :none, [:irfobs])
end

@testset "Test load_posterior_moments works" begin
    m <= Setting(:sampling_method, :MH)
    load_posterior_moments(m)
    moment_tables(m)
    moment_tables(m, use_mode = true)

    m <= Setting(:sampling_method, :SMC)
    m <= Setting(:n_parts, 20)
    m <= Setting(:n_Φ, 10)
    m <= Setting(:adaptive_tempering_target_smc, false)
    data = df_to_matrix(m, load_data(m))
    DSGE.smc2(m, data, verbose = :none, run_csminwel = false)
    load_posterior_moments(m)

end

@testset "Test meansbands_to_matrix works" begin
    meansbands_to_matrix(m, :full, :none, [:histobs])
end


fp = dirname(@__FILE__)
if isfile(joinpath(fp, "test.tex_forecast.tex"))
    rm(joinpath(fp, "test.tex_forecast.tex"))
end
if isfile(joinpath(fp, "test.tex_history.tex"))
    rm(joinpath(fp, "test.tex_history.tex"))
end
