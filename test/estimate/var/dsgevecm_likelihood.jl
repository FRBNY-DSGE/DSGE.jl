# See estimate/posterior.jl for a test of the wrapper for dsgevecm_likelihood that
# takes a DSGEVECM as an input.
using DSGE, FileIO, Test, ModelConstructors
fp = dirname(@__FILE__)
@testset "Likelihood function of DSGE-VECM p(Y | θ, λ)" begin
    matdata = load(joinpath(fp, "../../reference/test_dsgevecm_likelihood.jld2"))

    # note nobs in the Matlab code is actually T in the Julia code (n_observations)
    # and lambdaT in the Matlab code is λ * T in the Julia code.
    @test matdata["lnpy"] ≈ DSGE.dsgevecm_likelihood(matdata["YYYY"],
                                                     matdata["XXYY"], matdata["XXXX"],
                                                     matdata["yyyyd"], matdata["xxyyd"],
                                                     matdata["xxxxd"], Int(matdata["nobs"]),
                                                     matdata["lambdaT"] / matdata["nobs"],
                                                     Int(matdata["cointadd"]))
    @test matdata["lnpyinf"] ≈ DSGE.dsgevecm_likelihood(matdata["YYYY"],
                                                        matdata["XXYY"], matdata["XXXX"],
                                                        matdata["yyyyd"], matdata["xxyyd"],
                                                        matdata["xxxxd"], Int(matdata["nobs"]),
                                                        Inf,  Int(matdata["cointadd"]))
    @test_throws DomainError matdata["lnpyinf"] ≈ DSGE.dsgevecm_likelihood(matdata["YYYY"],
                                                                           matdata["XXYY"], matdata["XXXX"],
                                                                           matdata["yyyyd"], matdata["xxyyd"],
                                                                           matdata["xxxxd"], Int(matdata["nobs"]),
                                                                           -1., Int(matdata["cointadd"]))
end
