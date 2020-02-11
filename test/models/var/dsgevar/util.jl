fp = dirname(@__FILE__)
@testset "Draw stationary VAR coefficients and covariance matrix" begin
    matdata = load(joinpath(fp, "../../../reference/test_dsgevar_stationary_draws.jld2"))
    β_draw, Σ_draw = DSGE.draw_stationary_VAR(matdata["yyyyd"], matdata["xxyyd"],
                                         matdata["xxxxd"], Int(matdata["nobsc"]),
                                         Int(matdata["nvar"]), Int(matdata["nlags"]);
                                         test = true, test_Σ_draw = matdata["draw"],
                                         test_β_draw = vec(matdata["beta_draw"]))
    @test matdata["cct_sim"] ≈ β_draw
    @test matdata["sig_sim"] ≈ Σ_draw
end
