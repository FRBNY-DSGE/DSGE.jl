path = dirname(@__FILE__)

m = AnSchorfheide()

@testset "Test prepare forecast inputs" begin
    # Throw error if input_type=:subset but no subset_inds provided
    @test_throws ErrorException DSGE.prepare_forecast_inputs!(m, :subset, :none, [:histobs, :forecastobs])
    # Go back to normal policy rule
    m <= Setting(:alternative_policy, AltPolicy(:historical, eqcond, solve))
    # Test with df passed in
    if haskey(ENV, "FRED_API_KEY") || isfile(joinpath(homedir(), ".freddatarc"))
        df_in = load_data(m)
        global output_vars, df_out = DSGE.prepare_forecast_inputs!(m, :mode, :none, [:histobs, :forecastobs], df = df_in)
        @test output_vars == [:histobs, :forecastobs, :bddforecastobs]
        @test isapprox(df_to_matrix(m, df_in), df_to_matrix(m, df_out), atol = 0., nans = true)
        # Test without df passed in
        global output_vars, df_out = DSGE.prepare_forecast_inputs!(m, :mode, :none, [:histobs, :forecastobs])
        @test output_vars == [:histobs, :forecastobs, :bddforecastobs]
        @test isapprox(df_to_matrix(m, df_in), df_to_matrix(m, df_out), atol = 0., nans = true)
    else
        @warn "Skipping prepare forecast inputs test because FRED_API_KEY not present"
    end
end

@testset "Test load draws" begin
    global overrides = forecast_input_file_overrides(m)
    global overrides[:mode] = "$path/../reference/paramsmode_.h5"
    global overrides[:mode_draw_shocks] = "$path/../reference/paramsmode_.h5"

    # Init
    @test typeof(load_draws(m, :init, verbose = :none)) == Vector{Float64}
    @test typeof(load_draws(m, :init_draw_shocks, verbose = :none)) == Vector{Float64}
    @test typeof(load_draws(m, :prior, verbose = :none)) == Matrix{Float64}
    @test typeof(load_draws(m, :init_draw_shocks, 1:1, verbose = :none))==Vector{Vector{Float64}}
    @test typeof(load_draws(m, :prior, 1:1, verbose = :none))==Vector{Vector{Float64}}


    # MH
    m <= Setting(:sampling_method, :MH)
    global overrides[:full] = "$path/../reference/mhsave_.h5"
    # load_draws without blocks
    @test typeof(load_draws(m, :mode, verbose = :none)) == Vector{Float64}
    @test typeof(load_draws(m, :mode_draw_shocks, verbose = :none)) == Vector{Float64}
    @test typeof(load_draws(m, :full, verbose = :none)) == Matrix{Float64}
    @test typeof(load_draws(m, :subset, subset_inds = 1:10, verbose = :none)) == Matrix{Float64}

    # load_draws with blocks
    @test typeof(load_draws(m, :full, 1:1, verbose = :none))==Vector{Vector{Float64}}
    @test typeof(load_draws(m, :mode_draw_shocks, 1:1, verbose = :none))==Vector{Vector{Float64}}
       # can't block with empty/invalid block indices
    @test_throws ErrorException typeof(load_draws(m, :full, 1:0, verbose = :none))==Vector{Vector{Float64}}
       # can't block with only single draw
    @test_throws ErrorException load_draws(m, :mode, 1:1, verbose = :none)
    @test_throws ErrorException load_draws(m, :init, 1:1, verbose = :none)

    # SMC
    m <= Setting(:sampling_method, :SMC)
    global overrides[:full] = "$path/../reference/smcsave_.h5"
    # load draws without blocks
    @test typeof(load_draws(m, :mode, verbose = :none, use_highest_posterior_value = true)) == Vector{Float64}
    @test typeof(load_draws(m, :mode_draw_shocks, verbose = :none)) == Vector{Float64}
    @test typeof(load_draws(m, :full, verbose = :none)) == Matrix{Float64}
    @test typeof(load_draws(m, :subset, subset_inds = 1:10, verbose = :none, use_highest_posterior_value = true)) == Matrix{Float64}

    # load draws with blocks
    @test typeof(load_draws(m, :full, 1:1, verbose = :none))==Vector{Vector{Float64}}
    @test_throws ErrorException typeof(load_draws(m, :full, 1:0, verbose = :none))==Vector{Vector{Float64}}
    @test typeof(load_draws(m, :mode_draw_shocks, 1:1, verbose = :none))==Vector{Vector{Float64}}

    # If give nonsense sampler, should throw error
    m <= Setting(:sampling_method, :marco)
    @test_throws ErrorException load_draws(m, :mode, verbose = :none)

    # Test it can load the csminwel mode when using SMC
    delete!(overrides, :mode) # make sure no modal override
    m <= Setting(:sampling_method, :SMC)
    m <= Setting(:saveroot, joinpath(dirname(@__FILE__), "../reference"))
    m <= Setting(:data_vintage, "")
    @test typeof(load_draws(m, :mode, verbose = :none,
                            use_highest_posterior_value = true)) == Vector{Float64}
    @test typeof(load_draws(m, :mode, verbose = :none,
                            use_highest_posterior_value = false)) == Vector{Float64}
    m <= Setting(:fix, "true", true, "fix", "") # Add tail to estimation file name
    out_err = [false]
    try
        load_draws(m, :mode, verbose = :none, # Check no such smc_cloud file
                          use_highest_posterior_value = false)
    catch e
        if isa(e, SystemError) || isa(e, ArgumentError) # seems to be different depending on the Julia version
            out_err[1] = true
        else
            rethrow(e)
        end
    end
    @test out_err[1]
    try
        load_draws(m, :mode, verbose = :none, # Check no such smc_cloud file
                          use_highest_posterior_value = true)
    catch e
        if isa(e, SystemError) || isa(e, ArgumentError) # seems to be different depending on the Julia version
            out_err[1] = true
        else
            rethrow(e)
        end
    end
    @test out_err[1]
end
