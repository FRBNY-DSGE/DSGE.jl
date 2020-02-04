fp = dirname(@__FILE__)
data = load(joinpath(fp, "../../reference/var_irf_out.jld2"))
β = data["bhat_full"]
Σ = data["sigmahat_full"]
irfs1 = data["irfs_full_choleskyLR"]
irfs2 = data["irfs_full_maxBC"]
irfs3 = data["irfs_full_cholesky"]
out1 = impulse_responses(β, Σ, 1, 20;
                        method = :choleskyLR,
                        flip_shocks = true)
out2 = impulse_responses(β, Σ, 1, 20;
                        method = :maxBC,
                        flip_shocks = true)
out3 = impulse_responses(β, Σ, 1, 20;
                        method = :cholesky,
                        flip_shocks = true)

@testset "VAR impulse responses" begin
    @test @test_matrix_approx_eq irfs1 out1
    @test @test_matrix_approx_eq irfs2 out2
    @test @test_matrix_approx_eq irfs3 out3
end

nothing
