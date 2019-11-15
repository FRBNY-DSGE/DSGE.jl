using ModelConstructors, Nullables, SMC, Test, Distributed, Dates, DataFrames, OrderedCollections, FileIO, DataStructures, LinearAlgebra, StatsBase, Random, CSV, StateSpaceRoutines, HDF5, JLD2
import ModelConstructors: @test_matrix_approx_eq, @test_matrix_approx_eq_eps
@everywhere using DSGE, JLD2, Printf, LinearAlgebra, ModelConstructors, SMC

my_tests = [
            "core",
            "parameters",
            "models/an_schorfheide/an_schorfheide",
            "models/smets_wouters/smets_wouters",
            "models/smets_wouters_orig/smets_wouters_orig",
            "models/m990/m990",
            "models/m1002/m1002",
            "models/m1010/m1010",
            "models/poolmodel/poolmodel",
            "data/fred_data",
            "data/load_data",
            "data/load_data_poolmodel",
            "data/misc",
            "data/reverse_transform",
            "data/simulate_data",
            "data/transformations",
            "data/transform_data",
            "data/util",
            "solve/gensys",
            "solve/solve",
            "solve/solve_poolmodel",
            "estimate/filter",
            "estimate/cat",
            "estimate/posterior",
            "estimate/filter_poolmodel",
            "estimate/posterior_poolmodel",
            "estimate/estimate_bma",
            "estimate/hessian",
            "estimate/csminwel",
            "estimate/optimize",

            "estimate/metropolis_hastings",

            # "estimate/smc/smc",
            # "estimate/smc/helpers",
            # "estimate/smc/initialization",
            # "estimate/smc/util",
            # "estimate/smc/mutation",
            # "estimate/smc/resample",

            "forecast/drivers",
            "forecast/smooth",
            "forecast/forecast",
            "forecast/shock_decompositions",
            "forecast/impulse_responses",
            "forecast/io",
            "forecast/forecast_one",
            "forecast/util",
            "analysis/compute_meansbands",
            "analysis/df_to_table",
            "analysis/io",
            "analysis/meansbands",
            "analysis/moments",
            "analysis/util",
            "altpolicy/altpolicy",
            "altpolicy/taylor93",
            "altpolicy/taylor99",
            "scenarios/scenario",
            "scenarios/forecast",
            "scenarios/switching",
            "scenarios/drivers",
            "decomp/decompose_forecast",
            "plot/plot",
  	        "plot/util"
            #"packet/packet"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
