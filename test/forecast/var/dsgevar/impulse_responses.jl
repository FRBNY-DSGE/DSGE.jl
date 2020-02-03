using Test, DSGE, MAT

@testset "Impulse responses to identified impact matrix" begin
    matdata = matread("../../../reference/test_irfdsge.mat")
    irfout = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                               zeros(size(matdata["zz"], 1)), matdata["mmm"],
                               matdata["impact"], 1; accumulate = true,
                               cum_inds = 1)
    @test matdata["aairf"] ≈ irfout
end

@testset "Rotation of impulse responses to identified impact matrix" begin
    matdata1 = matread("../../../reference/test_irfdsge.mat")
    matdata2 = matread("../../../reference/test_rotations.mat")
    ŷ = DSGE.impulse_responses_rotation(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                                        zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                                        matdata1["impact"], Int(matdata2["k"]), matdata2["cct_sim"],
                                        matdata2["sig_sim"], matdata2["XXpred"],
                                        Int(matdata2["qahead"]); accumulate = true,
                               cum_inds = 1, test_shocks = matdata2["Shocks"])

    @test ŷ ≈ matdata2["yypred"]
end
