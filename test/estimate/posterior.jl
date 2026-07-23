using BenchmarkTools

path = dirname(@__FILE__)
@testset "Check prior overloading" begin
    custom_settings = [Setting(:date_forecast_start, quartertodate("2015-Q4"))]
    m = AnSchorfheide(custom_settings = custom_settings, testing = true)
    priordensity = exp(DSGE.prior(m))
    dsgevar = DSGE.DSGEVAR(m, collect(keys(m.exogenous_shocks)))
    @test prior(m) == prior(m.parameters)
    @test 0 <= priordensity <= 1 # ensure prior density is a density
    @test prior(dsgevar) == prior(dsgevar.dsge) == prior(m.parameters)
end

@testset "Check DSGE likelihood and posterior calculations" begin
    custom_settings = [Setting(:date_forecast_start, quartertodate("2015-Q4"))]
    m = AnSchorfheide(custom_settings = custom_settings, testing = true)

    file = "$path/../reference/posterior.jld2"
    data = Matrix{Float64}(load(file, "data")')
    lh_expected = load(file, "lh_expected")
    post_expected = load(file, "post_expected")

    lh = likelihood(m, data)
    @test lh_expected ≈ lh

    post = DSGE.posterior(m, data)
    @test post_expected ≈ post

    x = map(α->α.value, m.parameters)
    post_at_start = DSGE.posterior!(m, x, data)
    @test post_expected ≈ post_at_start

    # Ensure if we are not evaluating at start vector, then we do not get the reference
    # posterior
    global y = x .+ 0.01
    post_not_at_start = DSGE.posterior!(m, y, data)
    ϵ = 1.0
    @test abs(post_at_start - post_not_at_start) > ϵ

    # Check keywords
    m.parameters[2].value = 3.
    @test !isinf(likelihood(m, data; sampler = false))
    @test isinf(likelihood(m, data; sampler = true))
end

@testset "Check DSGEVAR likelihood and posterior calculations" begin
    dsge = Model1002("ss10")
    obs_i = [dsge.observables[:obs_nominalrate], dsge.observables[:obs_gdp],
             dsge.observables[:obs_gdpdeflator]]
    m = DSGE.DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)))
    DSGE.update!(m; observables = [:obs_nominalrate, :obs_gdp, :π_t], lags = 4, λ = 0.5)
    dsge_data =
        df_to_matrix(dsge,
                     CSV.read(joinpath(path, "../reference/test_dsgevar_likelihood_dsge_data.csv"), DataFrame))[obs_i, :]
    lh = likelihood(m, dsge_data)
    @test lh ≈ DSGE.dsgevar_likelihood(m, dsge_data) ≈
        load(joinpath(path, "../reference/test_dsgevar_likelihood_dsge.jld2"), "llh")

    post = DSGE.posterior(m, dsge_data)
    @test post ≈ lh + prior(m)

    x = map(α -> α.value, m.dsge.parameters)
    post_at_start = DSGE.posterior!(m, x, dsge_data)
    @test post ≈ post_at_start

    # Ensure if we are not evaluating at start vector, then we do not get the reference
    # posterior
    global y = x .+ 0.01
    post_not_at_start = DSGE.posterior!(m, y, dsge_data)
    ϵ = 1.0
    @test abs(post_at_start - post_not_at_start) > ϵ

    # Check keywords
    m.dsge.parameters[1].value = 3.
    @test !isinf(likelihood(m, dsge_data; sampler = false))
    @test isinf(likelihood(m, dsge_data; sampler = true))
end

################
# Benchmarking #
################
# Flip to true to run; off by default.
run_benchmarks = false

if run_benchmarks
    # DSGE (AnSchorfheide): prior / likelihood / posterior / posterior!.
    custom_settings = [Setting(:date_forecast_start, quartertodate("2015-Q4"))]
    m_dsge = AnSchorfheide(custom_settings = custom_settings, testing = true)
    data_dsge = Matrix{Float64}(load("$path/../reference/posterior.jld2", "data")')
    x_dsge = map(α -> α.value, m_dsge.parameters)

    b_prior = @benchmark prior($m_dsge)
    b_lh    = @benchmark likelihood($m_dsge, $data_dsge)
    b_post  = @benchmark DSGE.posterior($m_dsge, $data_dsge)
    b_postx = @benchmark DSGE.posterior!($m_dsge, $x_dsge, $data_dsge)

    # DSGE-VAR (Model1002 ss10): likelihood.
    dsge_var = Model1002("ss10")
    obs_i = [dsge_var.observables[:obs_nominalrate], dsge_var.observables[:obs_gdp],
             dsge_var.observables[:obs_gdpdeflator]]
    m_var = DSGE.DSGEVAR(dsge_var, collect(keys(dsge_var.exogenous_shocks)))
    DSGE.update!(m_var; observables = [:obs_nominalrate, :obs_gdp, :π_t], lags = 4, λ = 0.5)
    data_var = df_to_matrix(dsge_var,
                            CSV.read(joinpath(path, "../reference/test_dsgevar_likelihood_dsge_data.csv"),
                                     DataFrame))[obs_i, :]
    b_var_lh = @benchmark likelihood($m_var, $data_var)

    println("\n===== estimate/posterior benchmark results =====")
    for (name, b) in [("prior (dsge)       ", b_prior),
                      ("likelihood (dsge)  ", b_lh),
                      ("posterior (dsge)   ", b_post),
                      ("posterior! (dsge)  ", b_postx),
                      ("likelihood (dsgevar)", b_var_lh)]
        println(name, "  time:   ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end
