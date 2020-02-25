fp = dirname(@__FILE__)
@testset "DSGE impulse responses to a pre-specified impact matrix" begin
    matdata = load(joinpath(fp, "../../../reference/test_irfdsge.jld2"))
    irfout = impulse_responses(matdata["TTT"], matdata["RRR"], matdata["zz"],
                               zeros(size(matdata["zz"], 1)), matdata["mmm"],
                               matdata["impact"], 1; accumulate = true,
                               cum_inds = 1)
    nobs = size(matdata["zz"], 1)
    for i = 1:nobs
        @test vec(matdata["aairf"][:, nobs * (i - 1) + 1:nobs * i]) ≈
            vec(irfout[:, :, i])
    end
end

@testset "Impulse responses to structural shocks identified by a DSGE-VAR" begin
    matdata1 = load(joinpath(fp, "../../../reference/test_irfdsge.jld2"))
    matdata2 = load(joinpath(fp, "../../../reference/test_rotations.jld2"))
    ŷ = DSGE.impulse_responses(matdata1["TTT"], matdata1["RRR"], matdata1["zz"],
                               zeros(size(matdata1["zz"], 1)), matdata1["mmm"],
                               matdata1["impact"].^2, Int(matdata2["k"]), matdata2["cct_sim"],
                               matdata2["sig_sim"], vec(matdata2["XXpred"]),
                               Int(matdata2["qahead"]); accumulate = true,
                               cum_inds = 1, test_shocks =
                               convert(Matrix{Float64}, matdata2["Shocks"]'))

    @test ŷ ≈ matdata2["yypred"]
end
