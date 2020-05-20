fp = dirname(@__FILE__)

@testset "Impulse responses to structural shocks identified by a DSGE-VECM" begin
    matdata = load(joinpath(fp, "../../../reference/dsgevecm_lambda_irfs.jld2"))
    ŷ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                               vec(matdata["DD"]), matdata["MM"],
                               matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                               Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                               matdata["cointvec"], Int(matdata["qahead"]),
                               vec(matdata["XXpred"]);
                               test_shocks = matdata["Shocks"])

    @test ŷ ≈ matdata["yypred"]

    Random.seed!(1793) # drawing shocks, so this shouldn't be the same
    ŷ1 = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                vec(matdata["DD"]), matdata["MM"],
                                matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                matdata["cointvec"],
                                Int(matdata["qahead"]), vec(matdata["XXpred"]); draw_shocks = true)
    Random.seed!(1793) # test re-seeding works
    ŷ2 = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                vec(matdata["DD"]), matdata["MM"],
                                matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                matdata["cointvec"],
                                Int(matdata["qahead"]), vec(matdata["XXpred"]); draw_shocks = true)
    Random.seed!(1793) # flipping shocks should yield the same thing
    @info "The following warning is expected."
    ŷ3 = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                vec(matdata["DD"]), matdata["MM"],
                                matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                matdata["cointvec"],
                                Int(matdata["qahead"]), vec(matdata["XXpred"]); draw_shocks = true,
                                flip_shocks = true)

    @test !(ŷ1 ≈ matdata["yypred"])
    @test ŷ1 ≈ ŷ2
    @test ŷ1 ≈ ŷ3

    # Test drawing shocks
    Random.seed!(1793)
    ŷ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                               vec(matdata["DD"]), matdata["MM"],
                               matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                               Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                               matdata["cointvec"], 5, vec(matdata["XXpred"]);
                               draw_shocks = true)
    @test @test_matrix_approx_eq ŷ matdata["ypred_drawshocks"]

    # Now try passing in X̂ = zeros(...)
    test_shocks = matdata["Shocks"]
    test_shocks[:, 2:end] .= 0.
    ŷ1_fcast = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                      vec(matdata["DD"]), matdata["MM"],
                                      matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                      Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                      matdata["cointvec"], 5;
                                      test_shocks = test_shocks)

    # Compute forecast with no shock
    n_obs = Int(matdata["nvar"])
    p = 4
    baseline_ŷ = zeros(n_obs, 5) # baseline ŷ
    baseline_X̂ = zeros(length(matdata["XXpred"]))
    for t = 1:5
        out = vec(baseline_X̂' * matdata["cct_sim"])
        baseline_ŷ[:, t] = out
        addcoint = baseline_X̂[1:3] + matdata["cointvec"] * out
        baseline_X̂ .= vcat(addcoint, 1., out, baseline_X̂[3 + 2:Int(matdata["k"]) - n_obs])
    end
    ŷ1 = ŷ1_fcast  - baseline_ŷ

    # Now use deviations keyword
    ŷ2 = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                               vec(matdata["DD"]), matdata["MM"],
                               matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                               Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                               matdata["cointvec"], 5;
                               deviations = true,
                               test_shocks = test_shocks)

    # Replicate using VECM IRFs
    ZZ = matdata["ZZ"][1:n_obs, :]
    DD = matdata["DD"][1:n_obs]
    MM = matdata["MM"][1:n_obs, :]
    a0_m = convert(Matrix{Float64},
                   dropdims(impulse_responses(matdata["TTT"], matdata["RRR"], ZZ, DD, MM,
                                              sqrt.(matdata["QQ"]), 1,
                                              accumulate = false); dims = 2)')
    rotation, r_a = qr(a0_m)
    rotation = sign.(diag(r_a)) .* rotation'
    Σ_chol = cholesky(matdata["sig_sim"]).L * rotation

    expY = zeros(n_obs, p + 5)
    expY[:, p + 1] = Σ_chol * test_shocks[:, 1]
    addcoint = zeros(3)
    for t = 2:5
        addcoint += matdata["cointvec"] * expY[:, p + t - 1]
        xT = vcat(addcoint, vec(expY[:, p + t - 1:-1:t])) # stacks addcoint w/ lag yₜ₋₁, ..., yₜ₋ₚ
        expY[:, p + t] = vec(xT' * matdata["cct_sim"][vcat(1:3, 5:32), :])
    end
    expY = expY[:, p + 1:end]

    @test @test_matrix_approx_eq expY ŷ1
    @test @test_matrix_approx_eq expY ŷ2

    # Check IRFs to individual exogenous structural shocks
    ŷ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                               vec(matdata["DD"]), matdata["MM"],
                               matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                               Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                               matdata["cointvec"], 5, vec(matdata["XXpred"]))
    @test @test_matrix_approx_eq ŷ matdata["ypred_indivshocks"]

    # Now use deviations keyword
    ŷ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                               vec(matdata["DD"]), matdata["MM"],
                               matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                               Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                               matdata["cointvec"], 5, deviations = true)
    neg_ŷ = DSGE.impulse_responses(matdata["TTT"], matdata["RRR"], matdata["ZZ"],
                                   vec(matdata["DD"]), matdata["MM"],
                                   matdata["QQ"], Int(matdata["k"]), Int(matdata["nvar"]),
                                   Int(matdata["coint"]), matdata["cct_sim"], matdata["sig_sim"],
                                   matdata["cointvec"], 5, deviations = true,
                                   flip_shocks = true)
    @test @test_matrix_approx_eq neg_ŷ -ŷ

    # Replicate using VECM irfs
        ZZ = matdata["ZZ"][1:n_obs, :]
        DD = matdata["DD"][1:n_obs]
        MM = matdata["MM"][1:n_obs, :]
        a0_m = convert(Matrix{Float64},
                       dropdims(impulse_responses(matdata["TTT"], matdata["RRR"], ZZ, DD, MM,
                                                  sqrt.(matdata["QQ"]), 1,
                                                  accumulate = false); dims = 2)')
        rotation, r_a = qr(a0_m)
        rotation = sign.(diag(r_a)) .* rotation'
        Σ_chol = cholesky(matdata["sig_sim"]).L * rotation

    for i in 1:size(matdata["RRR"], 2)
        shocks = zeros(n_obs)
        shocks[i] = -sqrt(matdata["QQ"][i, i])
        expY = zeros(n_obs, p + 5)
        expY[:, p + 1] = Σ_chol * shocks
        addcoint = zeros(3)
        for t = 2:5
            addcoint += matdata["cointvec"] * expY[:, p + t - 1]
            xT = vcat(addcoint, vec(expY[:, p + t - 1:-1:t])) # stacks addcoint w/ lag yₜ₋₁, ..., yₜ₋ₚ
            expY[:, p + t] = vec(xT' * matdata["cct_sim"][vcat(1:3, 5:32), :])
        end
        expY = expY[:, p + 1:end]

        @test @test_matrix_approx_eq expY ŷ[:, :, i]
    end
end


nothing
