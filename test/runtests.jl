using Base.Test
using DSGE

include(joinpath(dirname(@__FILE__()),"util.jl"))

my_tests = [
            "core",
            "parameters",
            "models/an_schorfheide/an_schorfheide",
            "models/smets_wouters/smets_wouters",
            "models/m990/m990",
            "models/m1002/m1002",
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
            "estimate/eig",
            "estimate/metropolis_hastings",
            "forecast/smooth",
            "forecast/forecast",
            "forecast/shock_decompositions",
            "forecast/impulse_responses",
            "forecast/io",
            "forecast/forecast_one",
            "analysis/means_bands"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end