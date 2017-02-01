using Base.Test
using DSGE

include(joinpath(dirname(@__FILE__()),"util.jl"))

my_tests = [
            "core",
            "parameters",
            "models/m990/m990",
            "models/smets_wouters/smets_wouters",
            "data/misc",
            "data/load_data",
            "solve/gensys",
            "solve/solve",
            "estimate/kalman_filter",
            "estimate/kalman_filter_2part",
            "estimate/cat",
            "estimate/posterior",
            "estimate/hessizero",
            "estimate/hessian",
            "estimate/csminwel",
            "estimate/optimize",
            "estimate/eig",
            "estimate/metropolis_hastings",
            "forecast/smoothers",
            "forecast/filter",
            "forecast/filterandsmooth",
            "forecast/smooth",
            "forecast/forecast",
            "forecast/shock_decompositions",
            "forecast/impulse_responses",
            "forecast/forecast_one"
            ]

for test in my_tests
    test_file = string("$test.jl")
    @printf " * %s\n" test_file
    include(test_file)
end
