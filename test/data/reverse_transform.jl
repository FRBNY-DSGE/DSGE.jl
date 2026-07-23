using BenchmarkTools

# Load data to use for tests
path = dirname(@__FILE__)
fred = CSV.read("$path/../reference/fred_160812.csv", DataFrame)

# Specify vintage and dates
custom_settings = [Setting(:cond_id, 0),
                   Setting(:use_population_forecast, true),
                   Setting(:n_anticipated_shocks, 6),
                   Setting(:population_forecast, false)]
m = AnSchorfheide(testing = true, custom_settings = custom_settings)
m <= Setting(:saveroot, "$path/../reference/")
m <= Setting(:rate_expectations_source, :ois)

pop_forecast = DataFrame(CNP16OV = cumsum(fred[!,:CNP16OV][end] .* (ones(8) .* 1.03)))
pop_forecast[!,:date] = DSGE.get_quarter_ends(DSGE.next_quarter(fred[!,:date][end]), fred[!,:date][end] + Year(2))

mb_means = load("$path/../reference/means_bands_out.jld2")["exp_full_means"] # this is AnSchorfheide
m_dates = DSGE.quarter_range(DSGE.quartertodate("1960-Q1"), DSGE.quartertodate("2015-Q3"))
start_date = m_dates[1]
histobs_df = DataFrame(hcat(m_dates, mb_means[:histobs]'),
                       vcat(:date, collect(keys(m.observables))))
histpseudo_df = DataFrame(hcat(m_dates, mb_means[:histpseudo]'),
                          vcat(:date, collect(keys(m.pseudo_observables))))

revobs1 = Matrix{Float64}(reverse_transform(m, mb_means[:histobs], start_date,
                                            collect(keys(m.observables)), :obs, verbose = :none)[!,2:end])
revobs2 = Matrix{Float64}(reverse_transform(m, histobs_df, :obs, verbose = :none)[:,2:end])
revpseudo1 = Matrix{Float64}(reverse_transform(m, mb_means[:histpseudo], start_date,
                                               collect(keys(m.pseudo_observables)), :pseudo,
                                               verbose = :none)[:,2:end])
revpseudo2 = Matrix{Float64}(reverse_transform(m, histpseudo_df, :pseudo, verbose = :none)[!,2:end])
q4revobs1 = Matrix{Float64}(reverse_transform(m, mb_means[:histobs], start_date,
                                              collect(keys(m.observables)), :obs, verbose = :none,
                                              fourquarter = true)[4:end,2:end])
q4revobs2 = Matrix{Float64}(reverse_transform(m, histobs_df, :obs, verbose = :none,
                                              fourquarter = true)[4:end,2:end])
q4revpseudo1 = Matrix{Float64}(reverse_transform(m, mb_means[:histpseudo], start_date,
                                                 collect(keys(m.pseudo_observables)), :pseudo,
                                                 verbose = :none,
                                                 fourquarter = true)[4:end,2:end])
q4revpseudo2 = Matrix{Float64}(reverse_transform(m, histpseudo_df, :pseudo, verbose = :none,
                                                 fourquarter = true)[4:end,2:end])

exp_revobs, exp_revpseudo, exp_q4revobs, exp_q4revpseudo =
JLD2.jldopen("$path/../reference/reverse_transform_out.jld2", "r") do file
    read(file, "exp_revobs"), read(file, "exp_revpseudo"),
    read(file, "exp_q4revobs"), read(file, "exp_q4revpseudo")
end

@testset "Check reverse transform returns expected output" begin
    @test @test_matrix_approx_eq revobs1 exp_revobs
    @test @test_matrix_approx_eq revobs2 exp_revobs
    @test @test_matrix_approx_eq revpseudo1 exp_revpseudo
    @test @test_matrix_approx_eq revpseudo2 exp_revpseudo
    @test @test_matrix_approx_eq q4revobs1 exp_q4revobs
    @test @test_matrix_approx_eq q4revobs2 exp_q4revobs
    @test @test_matrix_approx_eq q4revpseudo1 exp_q4revpseudo
    @test @test_matrix_approx_eq q4revpseudo2 exp_q4revpseudo
end

################
# Benchmarking #
################
# Flip to true to run; off by default. None hit the FRED API.
run_benchmarks = false

if run_benchmarks
    obs_keys    = collect(keys(m.observables))
    pseudo_keys = collect(keys(m.pseudo_observables))

    # Matrix- and DataFrame-input variants (obs + pseudo), plus four-quarter.
    b_obs_mat = @benchmark reverse_transform($m, $(mb_means[:histobs]), $start_date,
                                             $obs_keys, :obs, verbose = :none)
    b_obs_df  = @benchmark reverse_transform($m, $histobs_df, :obs, verbose = :none)
    b_pseudo_mat = @benchmark reverse_transform($m, $(mb_means[:histpseudo]), $start_date,
                                                $pseudo_keys, :pseudo, verbose = :none)
    b_q4_obs  = @benchmark reverse_transform($m, $histobs_df, :obs, verbose = :none,
                                             fourquarter = true)

    println("\n===== data/reverse_transform benchmark results =====")
    for (name, b) in [("reverse_transform obs    (matrix)", b_obs_mat),
                      ("reverse_transform obs    (df)    ", b_obs_df),
                      ("reverse_transform pseudo (matrix)", b_pseudo_mat),
                      ("reverse_transform obs    (4q)    ", b_q4_obs)]
        println(name, "  time:   ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
