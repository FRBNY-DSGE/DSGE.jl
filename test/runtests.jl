using Test, Distributed, Dates, DataFrames, OrderedCollections, FileIO, DataStructures, LinearAlgebra, StatsBase, Random
@everywhere using DSGE, JLD2, Printf, LinearAlgebra

my_tests = [
            "core",
            "parameters",
            "models/an_schorfheide/an_schorfheide",
            "models/smets_wouters/smets_wouters",
            "models/m990/m990",
            "models/m1002/m1002",
            "models/m1010/m1010",
            "data/misc",
            "data/load_data",
            "solve/gensys",
            "solve/solve",
            "estimate/filter",
            "estimate/cat",
            "estimate/posterior",
            "estimate/hessizero",
            "estimate/hessian",
            "estimate/csminwel",
            "estimate/optimize",
            "estimate/metropolis_hastings",
            "data/misc",
            "data/load_data",
            "forecast/smooth",
            "forecast/forecast",
            "forecast/shock_decompositions",
            "forecast/impulse_responses",
            "forecast/io",
            "forecast/forecast_one",
            "analysis/compute_meansbands",
            "altpolicy/altpolicy",
            "scenarios/scenario",
            "scenarios/forecast",
            "scenarios/switching",
            "scenarios/drivers",
            "decomp/decompose_forecast",
            "plot/plot"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
