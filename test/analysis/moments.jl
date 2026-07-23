using DSGE, Test, FileIO, Random, ModelConstructors

# Set to true to regenerate the reference outputs in moments_poolmodel_outputs.jld2,
# then set back to false and rerun to confirm the test passes against them.
save_output = false

m = AnSchorfheide()
@testset "Test moments" begin
    @test DSGE.moments(m.parameters[findfirst(x -> x.key==:τ, m.parameters)]) == (2.0, 0.5)
    @test DSGE.moments(m.parameters[findfirst(x -> x.key==:e_y, m.parameters)]) == (0.1159846, 0.0)
    @test DSGE.moments(m.parameters[findfirst(x -> x.key==:σ_R, m.parameters)]) == (0.4, 4.0)
end

# This script currently just tests sample_λ, compute_Eλ
fp = dirname(@__FILE__)
out_path = "$(fp)/../reference/moments_poolmodel_outputs.jld2"

# --- Compute sample_λ outputs (dynamic pool) ---
Random.seed!(1793)
m = PoolModel("ss1")
θs = load("$(fp)/../reference/moments_poolmodel_inputs.jld2", "thetas")[1:2,:]
pred_dens = load("$(fp)/../reference/moments_poolmodel_inputs.jld2", "pred_dens")
λmat_noparallel = sample_λ(m, pred_dens, θs, 1)
λmat_parallel   = sample_λ(m, pred_dens, θs, 1; parallel = true)

# --- compute_Eλ (no reference outputs) ---
λvec = 0.5 * ones(2)
θchange = [.8 0. 1.] .* ones(2)
λhat_plush, λhat_t = compute_Eλ(m, 4, λvec, θchange)
λhat_plush_parallel, λhat_t_parallel = compute_Eλ(m, 4, λvec, θchange; parallel = true)

# --- Compute sample_λ outputs (static pool) ---
Random.seed!(1793)
sm = PoolModel("ss1"; weight_type = :static)
sm <= Setting(:saveroot, "$(fp)/../reference/")
sm <= Setting(:calculate_hessian, true)
sm <= Setting(:n_mh_simulations, 1)
sλmat_noparallel = sample_λ(sm, pred_dens, 1)
sλmat_parallel   = sample_λ(sm, pred_dens, 1; parallel = true)
try
    rm(joinpath(saveroot(sm), "output_data/poolmodel"); recursive = true, force = true)
catch err
    err isa Base.IOError || rethrow()
    @warn "Could not remove poolmodel output dir (likely a locked HDF5 file); restart Julia to release it." exception=err
end

# Optionally resave the reference outputs (merging into the existing file so any
# other stored keys are preserved).
if save_output
    ref = load(out_path)
    ref["lammat_noparallel"]      = λmat_noparallel
    ref["lammat_parallel"]        = λmat_parallel
    ref["slammat_noparallel_one"] = sλmat_noparallel
    ref["slammat_parallel_one"]   = sλmat_parallel
    save(out_path, ref)
end

# Load reference outputs (freshly saved if save_output, else the committed ones)
save_λmat_noparallel  = load(out_path, "lammat_noparallel")
save_λmat_parallel    = load(out_path, "lammat_parallel")
save_sλmat_noparallel = load(out_path, "slammat_noparallel_one")
save_sλmat_parallel   = load(out_path, "slammat_parallel_one")

@testset "Check sample_λ works correctly" begin
    @test @test_matrix_approx_eq save_λmat_noparallel λmat_noparallel
    @test @test_matrix_approx_eq save_λmat_parallel λmat_parallel
end

@testset "Check compute_Eλ works correctly" begin
    @test λhat_plush == 0.5
    @test λhat_plush_parallel == 0.5
    @test λhat_t == 0.5
    @test λhat_t_parallel == 0.5
end

@testset "Check sample_λ works for a static pool" begin
    @test sλmat_noparallel == save_sλmat_noparallel
    @test sλmat_parallel == save_sλmat_parallel
end
