using DSGE, Test, BenchmarkTools, Random

# Tests + benchmark for src/estimate/smc/util.jl.

module SMCUtilScaffold
    import DSGE: ParticleCloud, Cloud, VERBOSITY  # types for print signatures, VERBOSITY for their bodies
    using Random                                  # shuffle, in generate_free_blocks
    include(joinpath(@__DIR__, "..", "..", "..", "src", "estimate", "smc", "util.jl"))
end
smcutil = SMCUtilScaffold   # plain global alias (not const) so re-include never errors  

@testset "scalar_reduce" begin
    # n scalars per iteration -> n vectors collecting the i-th scalar across iterations
    r = smcutil.scalar_reduce([[1.0], [2.0]], [[3.0], [4.0]], [[5.0], [6.0]])
    @test r == [[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]]
end

@testset "vector_reduce" begin
    # n vectors per iteration -> n matrices, one column per iteration.
    # Inputs are the vector_reshape'd form (Vector{Matrix}), as in the SMC pipeline.
    i1 = smcutil.vector_reshape([1.0, 2.0])
    i2 = smcutil.vector_reshape([3.0, 4.0])
    i3 = smcutil.vector_reshape([5.0, 6.0])
    r = smcutil.vector_reduce(i1, i2, i3)
    @test length(r) == 1
    @test r[1] == [1.0 3.0 5.0; 2.0 4.0 6.0]
end

@testset "scalar_reshape" begin
    r = smcutil.scalar_reshape(1.0, 2.0, 3.0)
    @test r isa Vector{Vector{Float64}}
    @test r == [[1.0], [2.0], [3.0]]
    @test smcutil.scalar_reshape([1.0, 2.0], 3.0) == [[1.0, 2.0], [3.0]]   # vectors pass through
end

@testset "vector_reshape" begin
    r = smcutil.vector_reshape([1.0, 2.0, 3.0])
    @test r isa Vector{Matrix{Float64}}
    @test r[1] == reshape([1.0, 2.0, 3.0], 3, 1)
    r2 = smcutil.vector_reshape([1.0, 2.0], [3.0, 4.0, 5.0])
    @test size(r2[1]) == (2, 1)
    @test size(r2[2]) == (3, 1)
    @test smcutil.vector_reshape(5.0)[1] == reshape([5.0], 1, 1)   # scalar -> 1x1
end

@testset "generate_free_blocks" begin
    Random.seed!(47)
    blocks = smcutil.generate_free_blocks(10, 3)
    @test length(blocks) == 3
    @test sort(vcat(blocks...)) == collect(1:10)    # a partition of 1:10
    @test length.(blocks) == [4, 4, 2]              # ceil-div blocks, smaller last block

    one = smcutil.generate_free_blocks(5, 1)              # single block holds everything
    @test length(one) == 1
    @test sort(one[1]) == collect(1:5)

    singles = smcutil.generate_free_blocks(5, 5)          # one block per parameter
    @test length(singles) == 5
    @test all(length.(singles) .== 1)
    @test sort(vcat(singles...)) == collect(1:5)
end

@testset "generate_all_blocks" begin
    # map free-parameter block indices (1:n_free) back to full-parameter indices
    blocks_free    = [[1, 3], [2], [4, 5]]
    free_para_inds = [10, 20, 30, 40, 50]
    @test smcutil.generate_all_blocks(blocks_free, free_para_inds) == [[10, 30], [20], [40, 50]]
end

@testset "init/end_stage_print run without error (verbose=:low)" begin
    # verbose=:low skips the weighted-moment block, so a minimal cloud suffices.
    para_syms  = [:a, :b]
    mat_cloud  = DSGE.Cloud(zeros(4, 8), [0.5, 1.0], [100.0], 1, 2, 0, 0.5, 0.25, 12.0)
    part_cloud = DSGE.ParticleCloud(DSGE.Particle[], [0.5, 1.0], [100.0], 1, 2, 0, 0.5, 0.25, 12.0)
    # Suppress the print output. devnull isn't a redirectable stream on Julia 1.5,
    # so redirect into a temp file (a real IOStream) instead.
    mktemp() do _, io
        redirect_stdout(io) do
            @test smcutil.init_stage_print(part_cloud)            === nothing
            @test smcutil.init_stage_print(mat_cloud, para_syms)  === nothing
            @test smcutil.end_stage_print(part_cloud)             === nothing
            @test smcutil.end_stage_print(mat_cloud, para_syms)   === nothing
        end
    end
end

################
# Benchmarking #
################
run_benchmarks = false
if run_benchmarks
    blocks_free = [[1, 3], [2], [4, 5]]
    free_inds   = [10, 20, 30, 40, 50]

    b_gfb  = @benchmark smcutil.generate_free_blocks(50, 5)
    b_gab  = @benchmark smcutil.generate_all_blocks($blocks_free, $free_inds)
    b_sr   = @benchmark smcutil.scalar_reshape(1.0, 2.0, 3.0)
    b_vr   = @benchmark smcutil.vector_reshape([1.0, 2.0, 3.0])
    # reduce functions mutate their first arg, so build fresh inputs each sample
    b_sred = @benchmark smcutil.scalar_reduce(a, b, c) setup =
        (a = [[1.0], [2.0]]; b = [[3.0], [4.0]]; c = [[5.0], [6.0]])
    b_vred = @benchmark smcutil.vector_reduce(a, b, c) setup =
        (a = smcutil.vector_reshape([1.0, 2.0]); b = smcutil.vector_reshape([3.0, 4.0]); c = smcutil.vector_reshape([5.0, 6.0]))

    println("\n===== estimate/smc/util benchmark results =====")
    for (name, b) in [("generate_free_blocks", b_gfb),
                      ("generate_all_blocks", b_gab),
                      ("scalar_reshape", b_sr),
                      ("vector_reshape", b_vr),
                      ("scalar_reduce", b_sred),
                      ("vector_reduce", b_vred)]
        println(rpad(name, 22), " time: ", BenchmarkTools.prettytime(median(b).time),
                "   memory: ", BenchmarkTools.prettymemory(median(b).memory))
    end
end

nothing
