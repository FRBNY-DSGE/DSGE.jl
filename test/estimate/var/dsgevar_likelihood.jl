using DSGE, MAT, Test, CSV

@testset "Likelihood function of DSGE-VAR p(Y | θ, λ)" begin
    matdata = matread("../../reference/test_dsgevar_likelihood.mat")
    @test matdata["lnpy0"] ≈ DSGE.dsgevar_likelihood(matdata["YYYY"],
                                                matdata["XXYY"], matdata["XXXX"],
                                                matdata["yyyyd"], matdata["xxyyd"],
                                                matdata["xxxxd"], Int(matdata["nobs"]),
                                                matdata["lambda"] / matdata["nobs"],
                                                Int(matdata["nlags"]),
                                                Int(matdata["nvar"]))
    @test matdata["lnpyinf"] ≈ DSGE.dsgevar_likelihood(matdata["YYYY"],
                                                  matdata["XXYY"], matdata["XXXX"],
                                                  matdata["yyyyd"], matdata["xxyyd"],
                                                  matdata["xxxxd"], Int(matdata["nobs"]),
                                                  Inf, Int(matdata["nlags"]),
                                                  Int(matdata["nvar"]))


    dsge = Model1002("ss10")
    obs_i = [dsge.observables[:obs_nominalrate], dsge.observables[:obs_gdp], dsge.observables[:obs_gdpdeflator]]
    m = DSGE.DSGEVAR(dsge, collect(keys(dsge.exogenous_shocks)))
    DSGE.update!(m; observables = [:obs_nominalrate, :obs_gdp, :π_t], lags = 4)
    dsge_data = df_to_matrix(dsge, CSV.read("../../reference/test_dsgevar_likelihood_dsge_data.csv"))[obs_i, :]
    @test DSGE.dsgevar_likelihood(m, dsge_data; λ = .5) ≈ load("../../reference/test_dsgevar_likelihood_dsge.jld2", "llh")
end
