using FileIO
path = dirname(@__FILE__)

custom_settings = Dict{Symbol, Setting}(
    :date_forecast_start  => Setting(:date_forecast_start, quartertodate("2015-Q4")))
m = AnSchorfheide(custom_settings = custom_settings, testing = true)

file = "$path/../reference/posterior.jld2"
data = Matrix{Float64}(load(file, "data")')
lh_expected = load(file, "lh_expected")
post_expected = load(file, "post_expected")

@testset "Check likelihood and posterior calculations" begin
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
end
