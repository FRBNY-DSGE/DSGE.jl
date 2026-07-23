dsge = AnSchorfheide()

@testset "DSGEVAR measurement_error" begin
    m = DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)), "ss0")
    DSGE.update!(m; observables = [:obs_gdp, :obs_cpi, :obs_nominalrate], lags = 4)

    EE, MM = DSGE.measurement_error(m)

    @test EE == zeros(n_observables(m), n_observables(m))
    @test MM == zeros(n_observables(m), DSGE.n_shocks(m))
    @test size(EE) == (n_observables(m), n_observables(m))
    @test size(MM) == (n_observables(m), DSGE.n_shocks(m))
end

nothing
