path = dirname(@__FILE__)
@testset "Check prior overloading" begin
    custom_settings = Dict{Symbol, Setting}(
        :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
    m = AnSchorfheide(custom_settings = custom_settings, testing = true)
    priordensity = exp(DSGE.prior(m))
    dsgevar = DSGE.DSGEVAR(m, collect(keys(m.exogenous_shocks)))
    @test prior(m) == prior(m.parameters)
    @test 0 <= priordensity <= 1 # ensure prior density is a density
    @test prior(dsgevar) == prior(dsgevar.dsge) == prior(m.parameters)
end

@testset "Check DSGE likelihood and posterior calculations" begin
    custom_settings = Dict{Symbol, Setting}(
        :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
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
