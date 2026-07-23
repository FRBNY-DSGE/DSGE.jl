using DSGE, BenchmarkTools, Dates, DataFrames, OrderedCollections
path = dirname(@__FILE__)

# Construct a synthetic MeansBands in-memory to avoid stale JLD2 reference files.
# Band columns must use "LB"/"UB" patterns so which_density_bands recognises them.
let n = 20, vars = [:obs_gdp, :obs_cpi, :obs_nominalrate]
    dates = [Dates.lastdayofquarter(Date(2015, 10, 1) + Dates.Month(3*(i-1))) for i in 1:n]
    means_df = DataFrame(:date => dates, (v => Float64.(1:n) for v in vars)...)
    bands_dict = Dict{Symbol,DataFrame}(
        v => DataFrame(
            :date        => dates,
            Symbol("90.0% LB") => Float64.(1:n) .- 1.0,
            Symbol("90.0% UB") => Float64.(1:n) .+ 1.0,
        )
        for v in vars
    )
    metadata = Dict{Symbol,Any}(
        :product    => :forecast,
        :class      => :obs,
        :input_type => :full,
        :cond_type  => :none,
        :date_inds  => OrderedDict(d => i for (i, d) in enumerate(dates)),
        :indices    => OrderedDict(v => i for (i, v) in enumerate(vars)),
    )
    global mb_full = MeansBands(metadata, means_df, bands_dict)
end

means, bands = meansbands_to_matrix(mb_full)

nvars    = DSGE.n_vars_means(mb_full)
nperiods = DSGE.n_periods_means(mb_full)
nbands   = length(which_density_bands(mb_full))

@testset "meansbands_to_matrix" begin
    # Output shapes
    @test size(means) == (nvars, nperiods)
    @test size(bands) == (nbands, nvars, nperiods)

    # Output types
    @test eltype(means) == Float64
    @test eltype(bands) == Float64

    # Values match the input MeansBands means
    vars = DSGE.get_vars_means(mb_full)
    inds = mb_full.metadata[:indices]
    for v in vars
        @test means[inds[v], :] == mb_full.means[!, v]
    end
end

################
# Benchmarking #
################
run_benchmarks = false

if run_benchmarks
    b_m2m = @benchmark meansbands_to_matrix($mb_full)

    println("\n===== meansbands_to_matrix benchmark results =====")
    println(rpad("meansbands_to_matrix", 22), " time: ", rpad(BenchmarkTools.prettytime(median(b_m2m).time), 12),
            "memory: ", BenchmarkTools.prettymemory(median(b_m2m).memory))
end
