# See estimate/posterior.jl for a test of the wrapper for dsgevar_likelihood that
# takes a DSGEVAR as an input.
fp = dirname(@__FILE__)
@testset "Likelihood function of DSGE-VAR p(Y | θ, λ)" begin
    matdata = matread(joinpath(fp, "../../reference/test_dsgevar_likelihood.mat"))
    @test matdata["lnpy0"] ≈ DSGE.dsgevar_likelihood(matdata["YYYY"],
                                                matdata["XXYY"], matdata["XXXX"],
                                                matdata["yyyyd"], matdata["xxyyd"],
                                                matdata["xxxxd"], Int(matdata["nobs"]),
                                                matdata["lambda"] / matdata["nobs"])
    @test matdata["lnpyinf"] ≈ DSGE.dsgevar_likelihood(matdata["YYYY"],
                                                  matdata["XXYY"], matdata["XXXX"],
                                                  matdata["yyyyd"], matdata["xxyyd"],
                                                  matdata["xxxxd"], Int(matdata["nobs"]),
                                                  Inf)
    @test_throws DomainError matdata["lnpyinf"] ≈ DSGE.dsgevar_likelihood(matdata["YYYY"],
                                                  matdata["XXYY"], matdata["XXXX"],
                                                  matdata["yyyyd"], matdata["xxyyd"],
                                                  matdata["xxxxd"], Int(matdata["nobs"]),
                                                  -1.)
end
