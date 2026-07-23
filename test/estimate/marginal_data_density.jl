using DSGE, Test, BenchmarkTools, Random, JLD2, FileIO

# Tests + benchmark for src/estimate/marginal_data_density.jl.
#
# Unit-testable here: marginal_data_density(params, logpost, free_para_inds),
#   marginal_data_density_weighted(params, logpost, free_para_inds, cloud), tt2string.
#
# Integration here: marginal_data_density(m, data; :incremental_weights) drives the
#   default SMC path off a synthetic cloud .jld2 on disk (this is what online.jl
#   exercises). Not covered: the :harmonic_mean / :mh paths and
#   marginal_data_density_frontier (need real clouds or posterior! evaluations).

#-----------------------------------------------------------------
# core numeric: marginal_data_density(params, logpost, free_para_inds)
#-----------------------------------------------------------------
Random.seed!(47)
n_para, n_draws = 5, 400
params         = randn(n_para, n_draws)
logpost        = vec(-sum(params .^ 2, dims = 1) ./ 2 .- 3.0)   # finite, varied
free_para_inds = collect(1:n_para)

mdd = marginal_data_density(params, logpost, free_para_inds)

@testset "marginal_data_density (core numeric)" begin
    @test mdd isa Real
    @test isfinite(mdd)
end

#-----------------------------------------------------------------
# tt2string
#-----------------------------------------------------------------
@testset "tt2string" begin
    @test DSGE.tt2string(:new)   == "new"
    @test DSGE.tt2string(:old)   == "old"
    @test DSGE.tt2string(:whole) == "whole"
    @test DSGE.tt2string(:other) === nothing   # no else branch -> returns nothing
end

#-----------------------------------------------------------------
# weighted variant: marginal_data_density_weighted
# cloud is a particle matrix [para... loglh logprior old_loglh accept weight].
#-----------------------------------------------------------------
weights   = fill(1 / n_draws, n_draws)
cloud_mat = hcat(randn(n_draws, n_para), randn(n_draws), randn(n_draws),
                 randn(n_draws), fill(0.25, n_draws), weights)

wmdd = DSGE.marginal_data_density_weighted(params, logpost, free_para_inds, cloud_mat)

@testset "marginal_data_density_weighted" begin
    @test wmdd isa Real
    @test isfinite(wmdd)
end

#-----------------------------------------------------------------
# integration: marginal_data_density(m, data; :incremental_weights)
# Reads an SMC cloud .jld2 (keys cloud/w/W) from disk — the default path that
# online.jl exercises via marginal_data_density(m, data).
#-----------------------------------------------------------------
m          = AnSchorfheide(testing = true)
tmpdir     = mktempdir()
cloud_path = joinpath(tmpdir, "smc_cloud.jld2")

n_parts, n_stages = 100, 5
w_smc = rand(n_parts, n_stages) .+ 0.5     # positive -> finite column-sum logs
W_smc = rand(n_parts, n_stages) .+ 0.5
JLD2.jldopen(cloud_path, "w") do f
    f["cloud"] = zeros(2, 2)               # read but unused by the incremental path
    f["w"]     = w_smc
    f["W"]     = W_smc
end
dummy_data = randn(3, 20)                  # unused by the incremental path

@testset "marginal_data_density(m, data) :incremental_weights" begin
    out = marginal_data_density(m, dummy_data; estimation_method = :smc,
                                calculation_method = :incremental_weights,
                                smc_estimate_file = cloud_path)
    np       = sum(W_smc[:, 1])
    w_W      = w_smc[:, 2:end] .* W_smc[:, 1:end-1] ./ np
    expected = sum(log.(sum(w_W, dims = 1)))
    @test out ≈ expected
    @test isfinite(out)
end

@testset "marginal_data_density invalid method combo throws" begin
    @test_throws String marginal_data_density(m, dummy_data; estimation_method = :mh,
                                              calculation_method = :incremental_weights)
end

################
# Benchmarking #
################
run_benchmarks = false
if run_benchmarks
    b_mdd  = @benchmark marginal_data_density($params, $logpost, $free_para_inds)
    b_wmdd = @benchmark DSGE.marginal_data_density_weighted($params, $logpost,
                                                            $free_para_inds, $cloud_mat)
    b_tt   = @benchmark DSGE.tt2string(:whole)

    println("\n===== estimate/marginal_data_density benchmark results =====")
    for (name, b) in [("marginal_data_density (core)", b_mdd),
                      ("marginal_data_density_weighted", b_wmdd),
                      ("tt2string", b_tt)]
        println(rpad(name, 30), " time: ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
   
            end
end

nothing
