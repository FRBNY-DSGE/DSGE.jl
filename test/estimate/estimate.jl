using DSGE, Test, BenchmarkTools, HDF5, Statistics

# Unit tests + benchmark for the testable parts of src/estimate/estimate.jl:
#   - compute_parameter_covariance(path::String, method::Symbol; savepath)
#   - get_estimation_output_files(m)
#
# The rest of estimate.jl is integration-only and NOT unit-tested here:
#   - estimate(m, data::AbstractArray; ...) drives optimize!/hessian!/MH/SMC + H5
#     I/O and returns nothing; it's exercised end-to-end by the MH/SMC tests.
#   - estimate(m, df) / estimate(m) are thin wrappers around it (and currently
#     error with default args: their `old_data` default references a `data`
#     that isn't in scope — estimate.jl:83 and :113).
#   - compute_parameter_covariance(m; ...) just derives paths + delegates, so it
#     needs a model and a draws file on disk (light integration, not a unit).

#-----------------------------------------------------------------
# compute_parameter_covariance(path, method; savepath)
#-----------------------------------------------------------------
tmpdir = mktempdir()

@testset "compute_parameter_covariance from saved draws" begin
    for (method, prefix) in [(:MH, "mh"), (:SMC, "smc")]
        draws = randn(200, 6)                  # n_draws x n_params
        draws_path = joinpath(tmpdir, prefix * "save.h5")
        save_path  = joinpath(tmpdir, prefix * "_parameter_covariance.h5")
        h5open(draws_path, "w") do f
            f[prefix * "params"] = draws
        end

        compute_parameter_covariance(draws_path, method; savepath = save_path)

        @test isfile(save_path)
        saved = h5open(save_path, "r") do f
            read(f, prefix * "cov")
        end
        @test saved ≈ cov(draws)               # rows = observations
        @test size(saved) == (6, 6)
    end

    @testset "invalid method throws" begin
        draws_path = joinpath(tmpdir, "mhsave.h5")
        @test_throws String compute_parameter_covariance(draws_path, :NOT_A_METHOD)
    end

    @testset "missing draws file returns nothing" begin
        @test compute_parameter_covariance(joinpath(tmpdir, "absent.h5"), :MH) === nothing
    end
end

#-----------------------------------------------------------------
# get_estimation_output_files(m)
#-----------------------------------------------------------------
m = AnSchorfheide(testing = true)

@testset "get_estimation_output_files" begin
    files = get_estimation_output_files(m)

    @test files isa Dict{Symbol, String}
    @test Set(keys(files)) == Set([:paramsmode, :hessian, :mhsave, :paramsmean,
                                   :parameter_covariance, :priors,
                                   :prior_posterior_means, :moments])

    # Paths match the documented raw/work/tables locations.
    @test files[:paramsmode]            == rawpath(m, "estimate", "paramsmode.h5")
    @test files[:hessian]               == rawpath(m, "estimate", "hessian.h5")
    @test files[:mhsave]                == rawpath(m, "estimate", "mhsave.h5")
    @test files[:paramsmean]            == workpath(m, "estimate", "paramsmean.h5")
    @test files[:parameter_covariance]  == workpath(m, "estimate", "parameter_covariance.h5")
    @test files[:priors]                == tablespath(m, "estimate", "priors.tex")
    @test files[:prior_posterior_means] == tablespath(m, "estimate", "prior_posterior_means.tex")
    @test files[:moments]               == tablespath(m, "estimate", "moments.tex")
end

################
# Benchmarking #
################
# Flip to true to run; off by default. compute_parameter_covariance is disk I/O.
run_benchmarks = false

if run_benchmarks
    draws = randn(200, 6)
    draws_path = joinpath(tmpdir, "bench_mhsave.h5")
    save_path  = joinpath(tmpdir, "bench_cov.h5")
    h5open(draws_path, "w") do f
        f["mhparams"] = draws
    end

    b_cov   = @benchmark compute_parameter_covariance($draws_path, :MH; savepath = $save_path)
    b_files = @benchmark get_estimation_output_files($m)

    println("\n===== estimate/estimate benchmark results =====")
    for (name, b) in [("compute_parameter_covariance", b_cov),
                      ("get_estimation_output_files",  b_files)]
        println(rpad(name, 30), " time: ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
