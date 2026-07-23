using BenchmarkTools

file  = joinpath(dirname(@__FILE__), "smc", "reference", "smc_cloud_fix=true.jld2")
cloud = load(file, "cloud")

dsge_cloud = DSGE.Cloud(cloud)
smc_cloud  = SMC.Cloud(cloud)

# All three Cloud structs share these field names.
const CLOUD_FIELDS = (:particles, :tempering_schedule, :ESS, :stage_index,
                      :n_Φ, :resamples, :c, :accept, :total_sampling_time)

clouds_equal(a, b) = all(getfield(a, f) == getfield(b, f) for f in CLOUD_FIELDS)

N        = size(dsge_cloud.particles, 2)
n_params = DSGE.ind_para_end(N)
para_symbols = [Symbol("p", i) for i in 1:n_params]

@testset "backwards_compatibility: Cloud type converters" begin

    @testset "old_to_new_cloud(::DSGE.Cloud) -> SMC.Cloud" begin
        new = old_to_new_cloud(dsge_cloud)
        @test new isa SMC.Cloud
        @test clouds_equal(new, dsge_cloud)
    end

    @testset "SMC.Cloud(::DSGE.Cloud)" begin
        s = SMC.Cloud(dsge_cloud)
        @test s isa SMC.Cloud
        @test clouds_equal(s, dsge_cloud)
    end

    @testset "DSGE.Cloud(::SMC.Cloud)" begin
        d = DSGE.Cloud(smc_cloud)
        @test d isa DSGE.Cloud
        @test clouds_equal(d, smc_cloud)
    end

    @testset "SMC.Cloud(::SMC.Cloud) is the identity" begin
        @test SMC.Cloud(smc_cloud) === smc_cloud
    end

    @testset "round-trip DSGE.Cloud -> SMC.Cloud -> DSGE.Cloud" begin
        @test clouds_equal(DSGE.Cloud(SMC.Cloud(dsge_cloud)), dsge_cloud)
    end

    @testset "SMC.Cloud(::ParticleCloud) reconstructs the particle matrix" begin
        pc = ParticleCloud(dsge_cloud, para_symbols)
        @test pc isa ParticleCloud

        rebuilt = SMC.Cloud(pc)
        @test rebuilt isa SMC.Cloud
        # The particle matrix is rebuilt from per-particle fields; columns line up
        # exactly (params | loglh | logprior | old_loglh | accept | weight).
        @test rebuilt.particles ≈ dsge_cloud.particles
        # Scalar/vector metadata is carried through unchanged.
        for f in (:tempering_schedule, :ESS, :stage_index, :n_Φ,
                  :resamples, :c, :accept, :total_sampling_time)
            @test getfield(rebuilt, f) == getfield(dsge_cloud, f)
        end
    end
end

################
# Benchmarking #
################
# Flip to true to run; off by default. Inputs are local JLD2 data, no FRED API.
run_benchmarks = false

if run_benchmarks
    pc = ParticleCloud(dsge_cloud, para_symbols)

    b_oldnew   = @benchmark old_to_new_cloud($dsge_cloud)
    b_smc_dsge = @benchmark SMC.Cloud($dsge_cloud)
    b_dsge_smc = @benchmark DSGE.Cloud($smc_cloud)
    b_identity = @benchmark SMC.Cloud($smc_cloud)
    b_from_pc  = @benchmark SMC.Cloud($pc)

    println("\n===== estimate/smc/backwards_compatibility benchmark results =====")
    for (name, b) in [("old_to_new_cloud(::DSGE.Cloud)",  b_oldnew),
                      ("SMC.Cloud(::DSGE.Cloud)",          b_smc_dsge),
                      ("DSGE.Cloud(::SMC.Cloud)",          b_dsge_smc),
                      ("SMC.Cloud(::SMC.Cloud) identity",  b_identity),
                      ("SMC.Cloud(::ParticleCloud)",       b_from_pc)]
        println(rpad(name, 34), " time: ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
