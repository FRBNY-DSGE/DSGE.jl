using ModelConstructors, Nullables, SMC, Test, Distributed, Distributions
using Dates, DataFrames, OrderedCollections, FileIO, DataStructures, LinearAlgebra, SparseArrays
using StatsBase, Random, CSV, StateSpaceRoutines, HDF5, JLD2, MAT, Plots, Optim
import ModelConstructors: @test_matrix_approx_eq, @test_matrix_approx_eq_eps
@everywhere using DSGE, JLD2, Printf, LinearAlgebra, ModelConstructors, SMC

my_tests = [
    "parameters",
    "util",
    "defaults",
    "abstractdsgemodel",
    "abstractvarmodel",

    "altpolicy/altpolicy",
    "altpolicy/ait",
    "altpolicy/default_policy",
    "altpolicy/taylor93",

    "analysis/compute_meansbands",
    "analysis/create_q4q4_mb",
    "analysis/df_to_table",
    "analysis/io",
    "analysis/meansbands",
    "analysis/meansbands_to_matrix",
    "analysis/util",

    "data/fred_data",
    "data/reverse_transform",
    "data/simulate_data",
    "data/transformations",
    "data/transform_data",

    "decomp/decomposition_periods",

    "estimate/cat",
    "estimate/csminwel",
    "estimate/estimate",
    "estimate/filter",
    "estimate/hessian",
    "estimate/kalman",
    "estimate/optimize",
    "estimate/posterior",
    "estimate/resample",

    "estimate/smc/initialization",
    "estimate/smc/mutation",
    "estimate/smc/resample",
    "estimate/smc/particle",

    "estimate/var/dsgevar_likelihood",
    "estimate/var/dsgevecm_likelihood",

    "forecast/drivers",
    "forecast/smooth",
    "forecast/forecast",
    "forecast/shock_decompositions",
    "forecast/impulse_responses",
    "forecast/io",
    "forecast/forecast_one",
    "forecast/automatic_tempalt_zlb",
    "forecast/forecast_regime_switching",
    "forecast/time_varying_credibility",
    "forecast/multiple_altpol_imperfect_awareness",
    "forecast/m1002_ss62_forecast_test",
    "forecast/util",
    "forecast/var/impulse_responses",
    "forecast/var/dsgevar/impulse_responses",
    "forecast/var/dsgevecm/impulse_responses",

    "models/financial_frictions",
    "models/heterogeneous/bond_labor/bond_labor",
    "models/heterogeneous/het_dsge/het_dsge",
    "models/representative/an_schorfheide/an_schorfheide",
    "models/representative/m1002/m1002",
    "models/representative/m1010/m1010",

    "models/representative/smets_wouters/smets_wouters",
    "models/var/dsgevar/dsgevar",
    "models/var/dsgevecm/dsgevecm",

    "packet/packet",
    "plot/plot",

    "scenarios/scenario",
    "scenarios/forecast",
    "scenarios/switching",
    "scenarios/drivers",

    "solve/gensys",
    "solve/solve",
    "solve/gensys_uncertain_altpol",
    "solve/gensys2_uncertain_altpol_test1",
    "solve/gensys2_uncertain_altpol_test2",
    "solve/solve_poolmodel",

    "statespace/statespace",
]


failures = Tuple{String, Any}[]
for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    try
        include(test_file)
    catch err
        push!(failures, (test_file, err))
        @error "Test file failed" test_file exception = (err, catch_backtrace())
    end
end

# Full report across all test files (run to completion even if some throw).
println("\n", "="^70)
@printf "CI SUMMARY: %d of %d test files passed\n" (length(my_tests) - length(failures)) length(my_tests)
if !isempty(failures)
    println("Failed test files:")
    for (test_file, _) in failures
        println("  ✗ ", test_file)
    end
    error("$(length(failures)) test file(s) failed")
else
    println("All test files passed.")
end
