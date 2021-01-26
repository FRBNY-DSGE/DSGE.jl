using DSGE, FileIO, Random, Test, ModelConstructors
m = AnSchorfheide()
sys = compute_system(m)
path = dirname(@__FILE__)

@testset "Testing get_jstep " begin
    @test DSGE.get_jstep(m, 1)==1
    @test DSGE.get_jstep(m, 0)==5
end

@testset "Testing n_forecast_draws" begin
    @test DSGE.n_forecast_draws(m, :init) == 1
    @test DSGE.n_forecast_draws(m, :mode) == 1
    @test DSGE.n_forecast_draws(m, :mean) == 1

    m <= Setting(:sampling_method, :MH)
    global overrides = forecast_input_file_overrides(m)
    global overrides[:full] = "$path/../reference/mhsave_.h5"
    @test DSGE.n_forecast_draws(m, :full) == 100000

    m <= Setting(:sampling_method, :SMC)
    global overrides = forecast_input_file_overrides(m)
    global overrides[:full] = "$path/../reference/smcsave_.h5"
    @test DSGE.n_forecast_draws(m, :full) == 2000
    @test DSGE.n_forecast_draws(m, :prior) == 5000
    @test_throws ArgumentError DSGE.n_forecast_draws(m, :invalid_input)
end

@testset "Testing forecast_block_inds" begin
    m <= Setting(:jstep, 1)
    m <= Setting(:forecast_block_size, 5)
    block_inds, block_inds_thin = DSGE.forecast_block_inds(m, :full)
    @test block_inds[end] == 2000:5:2000
    @test block_inds_thin[end] == 400:400
    @test_throws ErrorException block_inds, block_inds_thin = DSGE.forecast_block_inds(m, :subset)
    block_inds, block_inds_thin = DSGE.forecast_block_inds(m, :subset, subset_inds = 1:100)
    @test block_inds[end] == 100:5:100
    @test block_inds_thin[end] == 20:20
    @test_throws ArgumentError DSGE.forecast_block_inds(m, :init)
    m <= Setting(:forecast_block_size, 5000)
    @test_throws ErrorException DSGE.forecast_block_inds(m, :full)
    @test_throws ErrorException DSGE.forecast_block_inds(m, :subset, subset_inds = 1:100)
    m <= Setting(:jstep, 5)
    m <= Setting(:forecast_block_size, 99)
    @test_throws ErrorException DSGE.forecast_block_inds(m, :full)
end

@testset "Testing add_requisite_output_vars" begin
    @test DSGE.add_requisite_output_vars([:histobs, :forecastobs, :shockdecobs, :shockdecpseudo]) == [:histobs, :forecastobs, :shockdecobs, :shockdecpseudo, :bddforecastobs, :dettrendobs, :dettrendpseudo, :trendobs, :trendpseudo]
    @test DSGE.add_requisite_output_vars([:histobs, :forecastobs]) == [:histobs, :forecastobs, :bddforecastobs]
    @test DSGE.add_requisite_output_vars([:shockdecobs]) == [:shockdecobs, :dettrendobs, :trendobs]
    @test_throws ErrorException DSGE.add_requisite_output_vars([:bad_input])
end

@testset "Testing remove_meansbands_only_output_vars" begin
    @test DSGE.remove_meansbands_only_output_vars([:histut, :hist4q, :forecastut, :forecast4q,
                                                   :bddforecastut, :bddforecast4q]) == []
    @test DSGE.remove_meansbands_only_output_vars([:histstates, :histut, :hist4q, :forecastut, :forecast4q,
                                                   :bddforecastut, :bddforecast4q]) == [:histstates]

end

@testset "Testing transplant history and forecast" begin
    @test isempty(DSGE.transplant_history(zeros(0, 0), 0))
    history = ones(10, 100)
    @test size(DSGE.transplant_history(history, 99)) == (10, 99)
    forecast = ones(10, 5)
    @test size(DSGE.transplant_forecast(history, forecast, 99)) == (10, 6)

    histstates = ones(8, 100)
    forecastobs = ones(3, 5)
    @test size(DSGE.transplant_forecast_observables(histstates, forecastobs, sys, 99)) == (3, 6)
end

@testset "Testing standardize_shocks" begin
    @test DSGE.standardize_shocks(ones(3, 3), sys[:QQ]) == ones(3,3) ./ sqrt.(diag(sys[:QQ]))
    QQs = [sys[:QQ], copy(sys[:QQ])]
    QQs[2][1, 1] = 1e3
    @test DSGE.standardize_shocks(ones(3, 3), QQs, [1:2, 3:3]) == hcat(ones(3, 2) ./ sqrt.(diag(sys[:QQ])),
                                                                     ones(3) ./ sqrt.(diag(QQs[2])))
end

if haskey(ENV, "FRED_API_KEY") || isfile(joinpath(homedir(),".freddatarc"))
    df = load_data(m)
    output_vars = [:histstates, :histobs, :forecaststates, :forecastobs]
    ndraws = 10
    paras = load_draws(m, :full, 1:ndraws, verbose = :none)

    forecast_outputs = [DSGE.forecast_one_draw(m, :full, :none, output_vars, p, df, verbose = :none) for p in paras]

    forecast_outputs = convert(Vector{Dict{Symbol, Array{Float64}}}, forecast_outputs)
    out = DSGE.assemble_block_outputs(forecast_outputs)
    @testset "Testing assemble_block_outputs" begin
        @test size(out[:histstates]) == (ndraws, n_states_augmented(m), n_mainsample_periods(m))
        @test size(out[:forecaststates]) == (ndraws, n_states_augmented(m), forecast_horizons(m))
        @test size(out[:forecastobs]) == (ndraws, n_observables(m), forecast_horizons(m))
    end

    @testset "Testing get_forecast_output_dims" begin
        @test DSGE.get_forecast_output_dims(m, :mode, :histobs)[1] == 1
        @test DSGE.get_forecast_output_dims(m, :mean, :histobs)[1] == 1
        @test DSGE.get_forecast_output_dims(m, :init, :histobs)[1] == 1
        m <= Setting(:forecast_block_size, 5)
        @test DSGE.get_forecast_output_dims(m, :full, :histobs)[1] == 400

        @test DSGE.get_forecast_output_dims(m, :mode, :histstates)[2] == n_states_augmented(m)
        @test DSGE.get_forecast_output_dims(m, :mode, :histobs)[2] == n_observables(m)
        @test DSGE.get_forecast_output_dims(m, :mode, :histpseudo)[2] == n_pseudo_observables(m)
        @test DSGE.get_forecast_output_dims(m, :mode, :histshocks)[2] == n_shocks_exogenous(m)

        @test DSGE.get_forecast_output_dims(m, :mode, :histstates)[3] == n_mainsample_periods(m)
        @test DSGE.get_forecast_output_dims(m, :mode, :forecaststates)[3] == forecast_horizons(m)

        @test DSGE.get_forecast_output_dims(m, :mode, :histshockdec)[3] == n_mainsample_periods(m)
        @test DSGE.get_forecast_output_dims(m, :mode, :irfobs)[3] == impulse_response_horizons(m)

        @test DSGE.get_forecast_output_dims(m, :mode, :trendobs) == (1, n_observables(m))

        @test DSGE.get_forecast_output_dims(m, :mode, :shockdecobs) == (1, n_observables(m), n_shockdec_periods(m), n_shocks_exogenous(m))
        @test DSGE.get_forecast_output_dims(m, :mode, :irfobs) == (1, n_observables(m), impulse_response_horizons(m), n_shocks_exogenous(m))
    end
else
    @warn "Skipping some tests because FRED_API_KEY not present"
end
